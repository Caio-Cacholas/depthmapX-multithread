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

#include "salalib/segmmodules/segmmetric.h"
#include "genlib/stringutils.h"
#include <thread>
#include <vector>
#include <atomic>
#include <cmath>
#include <algorithm>

bool SegmentMetric::run(Communicator* comm, ShapeGraph& map, bool) {

    AttributeTable& attributes = map.getAttributeTable();
    bool retvar = true;

    std::atomic<size_t> global_processed_rows(0);
    time_t atime = 0;

    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS,
            (m_sel_only ? map.getSelSet().size() : map.getConnections().size()));
    }
    int reccount = 0;

    size_t nShapes = map.getShapeCount();
    std::vector<int> axialrefs(nShapes);
    std::vector<float> seglengths(nShapes);

    int axial_col = attributes.getColumnIndex("Axial Line Ref");
    int length_col = attributes.getColumnIndex("Segment Length");

    // Descobre o número de threads
    // Descobre o número de threads
    int num_threads = num_threadsGlobal;
    if (num_threads <= 0) num_threads = 1;
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (max_threads == 0) max_threads = 1;
    if (num_threads > static_cast<int>(max_threads)) num_threads = static_cast<int>(max_threads);

    size_t chunk0 = (nShapes + num_threads - 1) / num_threads;

    std::vector<std::thread> threads0;

    // 1. Preencher axialrefs e seglengths em paralelo
    for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk0;
        size_t end = std::min(nShapes, (t + 1) * chunk0);
        threads0.emplace_back([&, start, end]() {
            for (size_t i = start; i < end; ++i) {
                AttributeRow& row = map.getAttributeRowFromShapeIndex(i);
                axialrefs[i] = row.getValue(axial_col);
                seglengths[i] = row.getValue(length_col);
            }
            });
    }
    for (auto& th : threads0) th.join();

    // 2. Calcular o máximo (redução)
    float maxseglength = *std::max_element(seglengths.begin(), seglengths.end());

    std::string prefix = "Metric ";
    std::string suffix = (m_radius != -1.0) ? dXstring::formatString(m_radius, " R%.f metric") : "";
    int maxbin = 512;

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

    // --- Multithread (C++11) ---
    struct Result {
        double totalmetdepth = 0.0, wtotal = 0.0, wtotaldepth = 0.0, totalsegdepth = 0.0;
        int total = 0;
        double meanDepth = 0.0, wmeanDepth = 0.0;
        double choice = 0.0, wchoice = 0.0;
        bool valid = false;
    };

    

    // Vetor para os resultados de cada thread
    std::vector<std::vector<Result>> results_per_thread(num_threads, std::vector<Result>(nShapes));

    auto process = [&](int start, int end, int thread_id) {
        time_t local_time = 0; // cada thread pode ter seu próprio timer
        for (int cursor = start; cursor < end; ++cursor) {
            AttributeRow& row = map.getAttributeRowFromShapeIndex(cursor);
            if (m_sel_only && !row.isSelected()) continue;

            std::vector<unsigned int> seen(nShapes, 0xffffffff);
            std::vector<TopoMetSegmentRef> audittrail(nShapes);
            std::vector<int> list[512]; // 512 bins!
            int bin = 0;
            list[bin].push_back(cursor);
            double rootseglength = seglengths[cursor];
            audittrail[cursor] = TopoMetSegmentRef(cursor, Connector::SEG_CONN_ALL, rootseglength * 0.5, -1);
            int open = 1;
            unsigned int segdepth = 0;
            double total = 0.0, wtotal = 0.0, wtotaldepth = 0.0, totalsegdepth = 0.0, totalmetdepth = 0.0;
            while (open != 0) {
                while (list[bin].size() == 0) {
                    bin++;
                    segdepth += 1;
                    if (bin == maxbin) {
                        bin = 0;
                    }
                }
                TopoMetSegmentRef& here = audittrail[list[bin].back()];
                list[bin].pop_back();
                open--;
                if (here.done) {
                    continue;
                }
                else {
                    here.done = true;
                }
                double len = seglengths[here.ref];
                totalsegdepth += segdepth;
                totalmetdepth += here.dist - len * 0.5; // preloaded with length ahead
                wtotal += len;
                wtotaldepth += len * (here.dist - len * 0.5);
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
                    if (!backsegs && iter == axline.m_forward_segconns.end()) {
                        break;
                    }

                    connected_cursor = iter->first.ref;

                    if (seen[connected_cursor] > segdepth && static_cast<size_t>(connected_cursor) != cursor) {
                        bool seenalready = (seen[connected_cursor] == 0xffffffff) ? false : true;
                        float length = seglengths[connected_cursor];
                        audittrail[connected_cursor] =
                            TopoMetSegmentRef(connected_cursor, here.dir, here.dist + length, here.ref);
                        seen[connected_cursor] = segdepth;
                        if (m_radius == -1 || here.dist + length < m_radius) {
                            open++;
                            list[(bin + int(floor(0.5 + 511 * length / maxseglength))) % 512].push_back(connected_cursor);
                        }
                        // Quick mod - TV
                        if (!m_sel_only && connected_cursor > int(cursor) && !seenalready) {
                            int subcur = connected_cursor;
                            while (subcur != -1) {
                                results_per_thread[thread_id][subcur].choice += 1;
                                results_per_thread[thread_id][subcur].wchoice += (rootseglength * length);
                                subcur = audittrail[subcur].previous;
                            }
                        }
                    }
                    iter++;
                }
            }

            // Preenche os resultados do thread
            if (total > 1) {
                Result& res = results_per_thread[thread_id][cursor];
                res.meanDepth = totalmetdepth / (total - 1);
                res.totalmetdepth = totalmetdepth;
                res.wmeanDepth = wtotaldepth / (wtotal - rootseglength);
                res.total = static_cast<int>(total);
                res.wtotal = wtotal;
                res.valid = true;
            }
            // === Barra de progresso por thread ===
            size_t processed = ++global_processed_rows;
            if (comm) {
                if (qtimer(local_time, 500)) { // Atualiza a cada 500ms por thread
                    if (comm->IsCancelled()) throw Communicator::CancelledException();
                    comm->CommPostMessage(Communicator::CURRENT_RECORD, processed);
                }
            }

        }
        };

    // Distribui as tarefas
    std::vector<std::thread> threads;
    int chunk = static_cast<int>(nShapes) / num_threads;
    for (int i = 0; i < num_threads; i++) {
        int start = i * chunk;
        int end = (i == num_threads - 1) ? static_cast<int>(nShapes) : start + chunk;
        threads.emplace_back(process, start, end, i);
    }

    for (auto& t : threads) t.join();

    // Merge dos resultados dos threads
    std::vector<Result> results(nShapes);
    for (int t = 0; t < num_threads; ++t) {
        for (size_t i = 0; i < nShapes; ++i) {
            results[i].choice += results_per_thread[t][i].choice;
            results[i].wchoice += results_per_thread[t][i].wchoice;
            results[i].meanDepth += results_per_thread[t][i].meanDepth;
            results[i].wmeanDepth += results_per_thread[t][i].wmeanDepth;
            results[i].totalmetdepth += results_per_thread[t][i].totalmetdepth;
            results[i].wtotal += results_per_thread[t][i].wtotal;
            results[i].total += results_per_thread[t][i].total;
            results[i].valid = results[i].valid || results_per_thread[t][i].valid;
        }
    }

    // Preenche os atributos
    for (size_t cursor = 0; cursor < nShapes; cursor++) {
        if (!results[cursor].valid) continue;
        AttributeRow& row = map.getAttributeRowFromShapeIndex(cursor);
        row.setValue(meandepthcol.c_str(), results[cursor].meanDepth);
        row.setValue(totaldcol.c_str(), results[cursor].totalmetdepth);
        row.setValue(wmeandepthcol.c_str(), results[cursor].wmeanDepth);
        row.setValue(totalcol.c_str(), results[cursor].total);
        row.setValue(wtotalcol.c_str(), results[cursor].wtotal);
    }

    if (!m_sel_only) {
        for (size_t cursor = 0; cursor < nShapes; cursor++) {
            AttributeRow& row = map.getAttributeRowFromShapeIndex(cursor);
            row.setValue(choicecol.c_str(), results[cursor].choice);
            row.setValue(wchoicecol.c_str(), results[cursor].wchoice);
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
