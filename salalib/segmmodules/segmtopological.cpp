// sala - a component of the depthmapX - spatial network analysis platform
// Copyright (C) 2000-2010, University College London, Alasdair Turner
// Copyright (C) 2011-2012, Tasos Varoudis
// Copyright (C) 2017-2018, Petros Koutsolampros

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <thread>
#include <vector>
#include <atomic>
#include <mutex>
//#include <QString>
#include <iostream>
//#include "depthmapX/settings.h"
//#include "depthmapX/settingsimpl.h"

#include "salalib/segmmodules/segmtopological.h"
#include "genlib/stringutils.h"

class TaskQueue {
public:
    TaskQueue(size_t max_tasks) : index(0), max(max_tasks) {}
    bool getNext(size_t& out) {
        size_t current = index.fetch_add(1);
        if (current < max) {
            out = current;
            return true;
        }
        return false;
    }
private:
    std::atomic<size_t> index;
    size_t max;
};

struct SegmentMetrics {
    bool valid = false;
    float meanDepth = 0.0f;
    float wmeanDepth = 0.0f;
    float totalDepth = 0.0f;
    int total = 0;
    float wtotal = 0.0f;

    float choice = 0.0f;
    float wchoice = 0.0f;

    std::mutex mutex;
};


void processSegmentTask(TaskQueue& queue,
                        ShapeGraph& map,
                        const std::vector<int>& axialrefs,
                        const std::vector<float>& seglengths,
                        std::vector<SegmentMetrics>& results,
                        bool m_sel_only,
                        double m_radius) {
    size_t cursor;
    while (queue.getNext(cursor)) {
        try {
            AttributeRow& row = map.getAttributeRowFromShapeIndex(cursor);
            if (m_sel_only && !row.isSelected()) continue;

            std::vector<unsigned int> seen(map.getShapeCount(), 0xffffffff);
            std::vector<TopoMetSegmentRef> audittrail(map.getShapeCount());
            std::vector<int> list[2];

            double rootseglength = seglengths[cursor];
            audittrail[cursor] = TopoMetSegmentRef(cursor, Connector::SEG_CONN_ALL, rootseglength * 0.5, -1);
            list[0].push_back(cursor);

            int open = 1, bin = 0, segdepth = 0;
            double total = 0.0, wtotal = 0.0, wtotaldepth = 0.0, totalsegdepth = 0.0;

            while (open != 0) {
                while (list[bin].empty()) {
                    bin = (bin + 1) % 2;
                    segdepth++;
                }

                TopoMetSegmentRef& here = audittrail[list[bin].back()];
                list[bin].pop_back(); open--;
                if (here.done) continue;
                here.done = true;

                double len = seglengths[here.ref];
                totalsegdepth += segdepth;
                wtotal += len;
                wtotaldepth += len * segdepth;
                total++;

                Connector& axline = map.getConnections().at(here.ref);
                int connected_cursor = -2;
                auto iter = axline.m_back_segconns.begin();
                bool backsegs = true;

                while (connected_cursor != -1) {
                    if (backsegs && iter == axline.m_back_segconns.end()) {
                        iter = axline.m_forward_segconns.begin(); backsegs = false;
                    }
                    if (!backsegs && iter == axline.m_forward_segconns.end()) break;

                    connected_cursor = iter->first.ref;
                    if (connected_cursor < 0 || static_cast<size_t>(connected_cursor) >= seglengths.size()) {
                        iter++;
                        continue;
                    }

                    if (seen[connected_cursor] > segdepth && size_t(connected_cursor) != cursor) {
                        float length = seglengths[connected_cursor];
                        int axialref = axialrefs[connected_cursor];

                        audittrail[connected_cursor] = TopoMetSegmentRef(connected_cursor, here.dir, here.dist + length, here.ref);
                        seen[connected_cursor] = segdepth;

                        if (m_radius == -1 || here.dist + length < m_radius) {
                            open++;
                            if (axialrefs[here.ref] == axialref) {
                                list[bin].push_back(connected_cursor);
                            } else {
                                list[(bin + 1) % 2].push_back(connected_cursor);
                                seen[connected_cursor] = segdepth + 1;
                            }
                        }

                        // Acumular choice apenas se não for seleção
                        if (!m_sel_only && connected_cursor > static_cast<int>(cursor)) {
                            int subcur = connected_cursor;
                            while (subcur != -1) {
                                std::lock_guard<std::mutex> lock(results[subcur].mutex);
                                results[subcur].choice += 1;
                                results[subcur].wchoice += rootseglength * length;
                                subcur = audittrail[subcur].previous;
                            }
                        }
                    }
                    iter++;
                }
            }

            if (total > 1) {
                SegmentMetrics& sm = results[cursor];
                sm.meanDepth = totalsegdepth / (total - 1);
                sm.totalDepth = totalsegdepth;
                sm.wmeanDepth = wtotaldepth / (wtotal - rootseglength);
                sm.total = static_cast<int>(total);
                sm.wtotal = wtotal;
                sm.valid = true;
            }

        } catch (...) {
            std::cerr << "Erro ao processar cursor " << cursor << "\n";
        }
    }
}

struct Result {
    float choice = 0;
    float wchoice = 0;
    float meanDepth = 0;
    float wmeanDepth = 0;
    float totalDepth = 0;
    float wtotal = 0;
    int total = 0;
    bool valid = false;
};

bool SegmentTopological::run(Communicator* comm, ShapeGraph& map, bool) {
    AttributeTable& attributes = map.getAttributeTable();
    bool retvar = true;
    std::atomic<size_t> global_processed_rows(0);
    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS,
            (m_sel_only ? map.getSelSet().size() : map.getConnections().size()));
    }

    size_t n = map.getShapeCount();
    std::vector<int> axialrefs(n);
    std::vector<float> seglengths(n);

    int axial_col = attributes.getColumnIndex("Axial Line Ref");
    int length_col = attributes.getColumnIndex("Segment Length");

    int num_threads = num_threadsGlobal;
    if (num_threads <= 0) num_threads = 1;

    unsigned int max_threads = std::thread::hardware_concurrency();
    if (max_threads == 0) max_threads = 1; // fallback

    if (num_threads > static_cast<int>(max_threads))
        num_threads = static_cast<int>(max_threads);

    size_t chunk0 = (n + num_threads - 1) / num_threads;
    std::vector<std::thread> threads0;

    // 1. Preencher vetores em paralelo
    for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk0;
        size_t end = std::min(n, (t + 1) * chunk0);
        threads0.emplace_back([&, start, end]() {
            for (size_t i = start; i < end; ++i) {
                AttributeRow& row = map.getAttributeRowFromShapeIndex(i);
                axialrefs[i] = row.getValue(axial_col);
                seglengths[i] = row.getValue(length_col);
            }
            });
    }
    for (auto& th : threads0) th.join();

    // 2. Calcular maxseglength (redução)
    float maxseglength = *std::max_element(seglengths.begin(), seglengths.end());

    std::string prefix = "Topological ";
    std::string suffix = (m_radius != -1.0) ? dXstring::formatString(m_radius, " R%.f metric") : "";
    std::string choicecol = prefix + "Choice" + suffix;
    std::string wchoicecol = prefix + "Choice [SLW]" + suffix;
    std::string meandepthcol = prefix + "Mean Depth" + suffix;
    std::string wmeandepthcol = prefix + "Mean Depth [SLW]" + suffix;
    std::string totaldcol = prefix + "Total Depth" + suffix;
    std::string totalcol = prefix + "Total Nodes" + suffix;
    std::string wtotalcol = prefix + "Total Length" + suffix;

    if (!m_sel_only) {
        attributes.insertOrResetColumn(choicecol.c_str());
        attributes.insertOrResetColumn(wchoicecol.c_str());
    }
    attributes.insertOrResetColumn(meandepthcol.c_str());
    attributes.insertOrResetColumn(wmeandepthcol.c_str());
    attributes.insertOrResetColumn(totaldcol.c_str());
    attributes.insertOrResetColumn(totalcol.c_str());
    attributes.insertOrResetColumn(wtotalcol.c_str());

    // --------- INÍCIO: OTIMIZAÇÃO PARA MULTITHREAD ---------
    

    size_t nShapes = map.getShapeCount();
    std::vector<std::vector<Result>> results_per_thread(num_threads, std::vector<Result>(nShapes));

    auto process = [&](int start, int end, int thread_id) {
        for (int cursor = start; cursor < end; ++cursor) {
            AttributeRow& row = map.getAttributeRowFromShapeIndex(cursor);
            if (m_sel_only && !row.isSelected()) continue;

            std::vector<unsigned int> seen(nShapes, 0xffffffff);
            std::vector<TopoMetSegmentRef> audittrail(nShapes);
            std::vector<int> list[2];
            int bin = 0;

            list[bin].push_back(cursor);
            double rootseglength = seglengths[cursor];
            audittrail[cursor] = TopoMetSegmentRef(cursor, Connector::SEG_CONN_ALL, rootseglength * 0.5, -1);

            int open = 1;
            unsigned int segdepth = 0;
            double total = 0.0, wtotal = 0.0, wtotaldepth = 0.0, totalsegdepth = 0.0;

            while (open != 0) {
                while (list[bin].empty()) {
                    bin = (bin + 1) % 2;
                    segdepth += 1;
                }

                TopoMetSegmentRef& here = audittrail[list[bin].back()];
                list[bin].pop_back();
                open--;

                if (here.done) continue;
                here.done = true;

                double len = seglengths[here.ref];
                totalsegdepth += segdepth;
                wtotal += len;
                wtotaldepth += len * segdepth;
                total += 1;

                Connector& axline = map.getConnections().at(here.ref);
                int connected_cursor = -2;

                auto iter = axline.m_back_segconns.begin();
                bool backsegs = true;

                while (connected_cursor != -1) {
                    if (backsegs && iter == axline.m_back_segconns.end()) {
                        iter = axline.m_forward_segconns.begin();
                        backsegs = false;
                    }
                    if (!backsegs && iter == axline.m_forward_segconns.end()) break;

                    connected_cursor = iter->first.ref;

                    if (connected_cursor != cursor && seen[connected_cursor] > segdepth) {
                        bool seenalready = (seen[connected_cursor] != 0xffffffff);
                        float length = seglengths[connected_cursor];
                        int axialref = axialrefs[connected_cursor];
                        audittrail[connected_cursor] = TopoMetSegmentRef(connected_cursor, here.dir, here.dist + length, here.ref);
                        seen[connected_cursor] = segdepth;

                        if (m_radius == -1 || here.dist + length < m_radius) {
                            open++;
                            if (axialrefs[here.ref] == axialref) {
                                list[bin].push_back(connected_cursor);
                            }
                            else {
                                list[(bin + 1) % 2].push_back(connected_cursor);
                                seen[connected_cursor] = segdepth + 1;
                            }
                        }

                        if (!m_sel_only && connected_cursor > cursor && !seenalready) {
                            int subcur = connected_cursor;
                            while (subcur != -1) {
                                // Sem mutex: cada thread escreve só no seu próprio array!
                                results_per_thread[thread_id][subcur].choice += 1;
                                results_per_thread[thread_id][subcur].wchoice += (rootseglength * length);
                                subcur = audittrail[subcur].previous;
                            }
                        }
                    }
                    ++iter;
                }
            }

            if (total > 1) {
                Result& r = results_per_thread[thread_id][cursor];
                r.meanDepth = totalsegdepth / (total - 1);
                r.totalDepth = totalsegdepth;
                r.wmeanDepth = wtotaldepth / (wtotal - rootseglength);
                r.total = static_cast<int>(total);
                r.wtotal = wtotal;
                r.valid = true;
            }
            size_t processed = ++global_processed_rows;  // incrementa e retorna o valor atual
            if (comm && qtimer(atime, 500)) {
                if (comm->IsCancelled()) throw Communicator::CancelledException();
                comm->CommPostMessage(Communicator::CURRENT_RECORD, processed);
            }
        }
    };

    std::vector<std::thread> threads;
    int chunk = static_cast<int>(nShapes) / num_threads;
    for (int i = 0; i < num_threads; i++) {
        int start = i * chunk;
        int end = (i == num_threads - 1) ? static_cast<int>(nShapes) : start + chunk;
        threads.emplace_back(process, start, end, i);
    }

    for (auto& t : threads) t.join();

    // ------ FAZER O MERGE DE TODOS OS RESULTADOS ------
    std::vector<Result> results(nShapes);
    for (int t = 0; t < num_threads; ++t) {
        for (size_t i = 0; i < nShapes; ++i) {
            results[i].choice += results_per_thread[t][i].choice;
            results[i].wchoice += results_per_thread[t][i].wchoice;
            results[i].meanDepth += results_per_thread[t][i].meanDepth;
            results[i].wmeanDepth += results_per_thread[t][i].wmeanDepth;
            results[i].totalDepth += results_per_thread[t][i].totalDepth;
            results[i].wtotal += results_per_thread[t][i].wtotal;
            results[i].total += results_per_thread[t][i].total;
            results[i].valid = results[i].valid || results_per_thread[t][i].valid;
        }
    }
    // --------- FIM: OTIMIZAÇÃO PARA MULTITHREAD ---------

    // Atualiza os atributos normalmente usando results[]
    for (size_t cursor = 0; cursor < nShapes; cursor++) {
        if (!results[cursor].valid) continue;
        AttributeRow& row = map.getAttributeRowFromShapeIndex(cursor);
        row.setValue(attributes.getColumnIndex(meandepthcol.c_str()), results[cursor].meanDepth);
        row.setValue(attributes.getColumnIndex(totaldcol.c_str()), results[cursor].totalDepth);
        row.setValue(attributes.getColumnIndex(wmeandepthcol.c_str()), results[cursor].wmeanDepth);
        row.setValue(attributes.getColumnIndex(totalcol.c_str()), static_cast<float>(results[cursor].total));
        row.setValue(attributes.getColumnIndex(wtotalcol.c_str()), results[cursor].wtotal);
    }

    if (!m_sel_only) {
        for (size_t cursor = 0; cursor < nShapes; cursor++) {
            AttributeRow& row = map.getAttributeRowFromShapeIndex(cursor);
            row.setValue(attributes.getColumnIndex(choicecol.c_str()), results[cursor].choice);
            row.setValue(attributes.getColumnIndex(wchoicecol.c_str()), results[cursor].wchoice);
        }
    }

    if (!m_sel_only) {
        map.setDisplayedAttribute(attributes.getColumnIndex(choicecol.c_str()));
    }
    else {
        map.setDisplayedAttribute(attributes.getColumnIndex(meandepthcol.c_str()));
    }

    return retvar;
}

