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

#include "salalib/segmmodules/segmtulip.h"

#include "genlib/stringutils.h"

//#include <QString>

//#include "depthmapX/settings.h"
//#include "depthmapX/settingsimpl.h"

#include <thread>
#include <mutex>
#include <atomic>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <array>

struct AuditTrailCell {
    int depth = -1;
    SegmentRef previous;
    bool leaf = true;
    void clearLine() { depth = -1; previous = SegmentRef(); leaf = true; }
};



bool SegmentTulip::run(Communicator* comm, ShapeGraph& map, bool) {
    std::atomic<size_t> global_processed_rows(0);

    if (map.getMapType() != ShapeMap::SEGMENTMAP) {
        return false;
    }

    int weighting_col2 = m_weighted_measure_col2;
    int routeweight_col = m_routeweight_col;
    bool interactive = m_interactive;
    AttributeTable& attributes = map.getAttributeTable();
    //int processed_rows = 0;
    time_t atime = 0;

    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS,
            (m_sel_only ? map.getSelSet().size() : map.getConnections().size()));
    }

    // Processamento dos raios e variáveis auxiliares
    bool radius_n = false;
    std::vector<double> radius_unconverted;
    for (int radius : m_radius_set) {
        if (radius == -1.0) radius_n = true;
        else radius_unconverted.push_back(radius);
    }
    if (radius_n) radius_unconverted.push_back(-1.0);

    std::vector<float> weights, routeweights;
    std::string weighting_col_text;
    int tulip_bins = m_tulip_bins;

    // Prepara os vetores de peso
    if (m_weighted_measure_col != -1) {
        weighting_col_text = attributes.getColumnName(m_weighted_measure_col);
        for (size_t i = 0; i < map.getConnections().size(); i++)
            weights.push_back(map.getAttributeRowFromShapeIndex(i).getValue(m_weighted_measure_col));
    }
    else {
        weights.assign(map.getConnections().size(), 1.0f);
    }
    std::string routeweight_col_text;
    if (routeweight_col != -1) {
        double max_value = attributes.getColumn(routeweight_col).getStats().max;
        routeweight_col_text = attributes.getColumnName(routeweight_col);
        for (size_t i = 0; i < map.getConnections().size(); i++)
            routeweights.push_back(1.0 - (map.getAttributeRowFromShapeIndex(i).getValue(routeweight_col) / max_value));
    }
    else {
        routeweights.assign(map.getConnections().size(), 1.0f);
    }

    // EFEF*
    std::vector<float> weights2;
    std::string weighting_col_text2;
    if (weighting_col2 != -1) {
        weighting_col_text2 = attributes.getColumnName(weighting_col2);
        for (size_t i = 0; i < map.getConnections().size(); i++)
            weights2.push_back(map.getAttributeRowFromShapeIndex(i).getValue(weighting_col2));
    }
    else {
        weights2.assign(map.getConnections().size(), 1.0f);
    }
    //*EFEF

    std::string tulip_text = std::string("T") + dXstring::formatString(tulip_bins, "%d");

    // Insere e obtém índices das colunas de atributo necessárias
    size_t r;
    for (r = 0; r < radius_unconverted.size(); r++) {
        std::string radius_text = makeRadiusText(m_radius_type, radius_unconverted[r]);
        if (m_choice) {
            if (routeweight_col != -1) {
                std::string choice_col_text = tulip_text + " Choice [Route weight by " + routeweight_col_text + "]" + radius_text;
                attributes.insertOrResetColumn(choice_col_text.c_str());
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text = tulip_text + " Choice [[Route weight by " + routeweight_col_text + "][" + weighting_col_text + " Wgt]]" + radius_text;
                    attributes.insertOrResetColumn(w_choice_col_text.c_str());
                }
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [[Route weight by " + routeweight_col_text + "][" + weighting_col_text + "-" + weighting_col_text2 + " Wgt]]" + radius_text;
                    attributes.insertOrResetColumn(w_choice_col_text2.c_str());
                }
            }
            else {
                std::string choice_col_text = tulip_text + " Choice" + radius_text;
                attributes.insertOrResetColumn(choice_col_text.c_str());
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text = tulip_text + " Choice [" + weighting_col_text + " Wgt]" + radius_text;
                    attributes.insertOrResetColumn(w_choice_col_text.c_str());
                }
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [" + weighting_col_text + "-" + weighting_col_text2 + " Wgt]" + radius_text;
                    attributes.insertOrResetColumn(w_choice_col_text2.c_str());
                }
            }
        }
        // Após inserir as colunas de Choice acima, adicione as de Integration, Node Count, Total Depth, etc:
        if (routeweight_col != -1) {
            std::string integ_col_text = tulip_text + " Integration [Route weight by " + routeweight_col_text + "]" + radius_text;
            std::string w_integ_col_text = tulip_text + " Integration [[Route weight by " + routeweight_col_text +
                "][" + weighting_col_text + " Wgt]]" + radius_text;

            std::string count_col_text = tulip_text + " Node Count [Route weight by " + routeweight_col_text + "]" + radius_text;
            std::string td_col_text = tulip_text + " Total Depth [Route weight by " + routeweight_col_text + "]" + radius_text;
            std::string w_td_text = tulip_text + " Total Depth [[Route weight by " + routeweight_col_text + "][" +
                weighting_col_text + " Wgt]]" + radius_text;
            std::string total_weight_text = tulip_text + " Total " + weighting_col_text + " [Route weight by " +
                routeweight_col_text + "]" + radius_text;

            attributes.insertOrResetColumn(integ_col_text.c_str());
            attributes.insertOrResetColumn(count_col_text.c_str());
            attributes.insertOrResetColumn(td_col_text.c_str());
            if (m_weighted_measure_col != -1) {
                attributes.insertOrResetColumn(w_integ_col_text.c_str());
                attributes.insertOrResetColumn(w_td_text.c_str());
                attributes.insertOrResetColumn(total_weight_text.c_str());
            }
        }
        else {
            std::string integ_col_text = tulip_text + " Integration" + radius_text;
            std::string w_integ_col_text = tulip_text + " Integration [" + weighting_col_text + " Wgt]" + radius_text;

            std::string count_col_text = tulip_text + " Node Count" + radius_text;
            std::string td_col_text = tulip_text + " Total Depth" + radius_text;
            std::string w_td_text = tulip_text + " Total Depth [" + weighting_col_text + " Wgt]" + radius_text;
            std::string total_weight_text = tulip_text + " Total " + weighting_col_text + radius_text;

            attributes.insertOrResetColumn(integ_col_text.c_str());
            attributes.insertOrResetColumn(count_col_text.c_str());
            attributes.insertOrResetColumn(td_col_text.c_str());
            if (m_weighted_measure_col != -1) {
                attributes.insertOrResetColumn(w_integ_col_text.c_str());
                attributes.insertOrResetColumn(w_td_text.c_str());
                attributes.insertOrResetColumn(total_weight_text.c_str());
            }
        }
    }
    // (Obtém os índices das colunas para cada tipo igual ao original)

    std::vector<int> choice_col, w_choice_col, w_choice_col2, count_col, integ_col, w_integ_col, td_col, w_td_col, total_weight_col;

    // Para cada radius, obtenha os índices das colunas, igual ao seu código:
    for (size_t r = 0; r < radius_unconverted.size(); r++) {
        std::string radius_text = makeRadiusText(m_radius_type, radius_unconverted[r]);
        if (m_choice) {
            if (routeweight_col != -1) {
                std::string choice_col_text = tulip_text + " Choice [Route weight by " + routeweight_col_text + "]" + radius_text;
                choice_col.push_back(attributes.getColumnIndex(choice_col_text.c_str()));
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text = tulip_text + " Choice [[Route weight by " + routeweight_col_text +
                        "][" + weighting_col_text + " Wgt]]" + radius_text;
                    w_choice_col.push_back(attributes.getColumnIndex(w_choice_col_text.c_str()));
                }
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [[Route weight by " + routeweight_col_text +
                        "][" + weighting_col_text + "-" + weighting_col_text2 + " Wgt]]" + radius_text;
                    w_choice_col2.push_back(attributes.getColumnIndex(w_choice_col_text2.c_str()));
                }
            }
            else {
                std::string choice_col_text = tulip_text + " Choice" + radius_text;
                choice_col.push_back(attributes.getColumnIndex(choice_col_text.c_str()));
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text = tulip_text + " Choice [" + weighting_col_text + " Wgt]" + radius_text;
                    w_choice_col.push_back(attributes.getColumnIndex(w_choice_col_text.c_str()));
                }
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [" + weighting_col_text + "-" +
                        weighting_col_text2 + " Wgt]" + radius_text;
                    w_choice_col2.push_back(attributes.getColumnIndex(w_choice_col_text2.c_str()));
                }
            }
        }
        if (routeweight_col != -1) {
            std::string integ_col_text = tulip_text + " Integration [Route weight by " + routeweight_col_text + "]" + radius_text;
            std::string w_integ_col_text = tulip_text + " Integration [[Route weight by " + routeweight_col_text +
                "][" + weighting_col_text + " Wgt]]" + radius_text;

            std::string count_col_text = tulip_text + " Node Count [Route weight by " + routeweight_col_text + "]" + radius_text;
            std::string td_col_text = tulip_text + " Total Depth [Route weight by " + routeweight_col_text + "]" + radius_text;
            std::string w_td_text = tulip_text + " Total Depth [[Route weight by " + routeweight_col_text + "][" +
                weighting_col_text + " Wgt]]" + radius_text;
            std::string total_weight_col_text = tulip_text + " Total " + weighting_col_text + " [Route weight by " +
                routeweight_col_text + "]" + radius_text;

            integ_col.push_back(attributes.getColumnIndex(integ_col_text.c_str()));
            count_col.push_back(attributes.getColumnIndex(count_col_text.c_str()));
            td_col.push_back(attributes.getColumnIndex(td_col_text.c_str()));
            if (m_weighted_measure_col != -1) {
                w_integ_col.push_back(attributes.getColumnIndex(w_integ_col_text.c_str()));
                w_td_col.push_back(attributes.getColumnIndex(w_td_text.c_str()));
                total_weight_col.push_back(attributes.getColumnIndex(total_weight_col_text.c_str()));
            }
        }
        else {
            std::string integ_col_text = tulip_text + " Integration" + radius_text;
            std::string w_integ_col_text = tulip_text + " Integration [" + weighting_col_text + " Wgt]" + radius_text;

            std::string count_col_text = tulip_text + " Node Count" + radius_text;
            std::string td_col_text = tulip_text + " Total Depth" + radius_text;
            std::string w_td_text = tulip_text + " Total Depth [" + weighting_col_text + " Wgt]" + radius_text;
            std::string total_weight_col_text = tulip_text + " Total " + weighting_col_text + radius_text;

            integ_col.push_back(attributes.getColumnIndex(integ_col_text.c_str()));
            count_col.push_back(attributes.getColumnIndex(count_col_text.c_str()));
            td_col.push_back(attributes.getColumnIndex(td_col_text.c_str()));
            if (m_weighted_measure_col != -1) {
                w_integ_col.push_back(attributes.getColumnIndex(w_integ_col_text.c_str()));
                w_td_col.push_back(attributes.getColumnIndex(w_td_text.c_str()));
                total_weight_col.push_back(attributes.getColumnIndex(total_weight_col_text.c_str()));
            }
        }
    }

    tulip_bins /= 2;
    tulip_bins += 1;

    std::vector<double> radius;
    for (r = 0; r < radius_unconverted.size(); r++) {
        if (m_radius_type == Options::RADIUS_ANGULAR && radius_unconverted[r] != -1)
            radius.push_back(floor(radius_unconverted[r] * tulip_bins * 0.5));
        else
            radius.push_back(radius_unconverted[r]);
    }

    int num_threads = num_threadsGlobal;
    if (num_threads <= 0) num_threads = 1;
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (max_threads == 0) max_threads = 1;
    if (num_threads > static_cast<int>(max_threads)) num_threads = static_cast<int>(max_threads);

    // --- Pré-cálculo do comprimento dos segmentos (threaded) ---
    int length_col = attributes.getColumnIndex("Segment Length");
    size_t n = map.getConnections().size();
    std::vector<float> lengths(n, 1.0f);
    if (length_col != -1) {
        size_t chunk = (n + num_threads - 1) / num_threads;
        std::vector<std::thread> threads;
        for (unsigned int t = 0; t < num_threads; ++t) {
            size_t start = t * chunk;
            size_t end = std::min(n, (t + 1) * chunk);
            threads.emplace_back([&, start, end]() {
                for (size_t i = start; i < end; ++i)
                    lengths[i] = map.getAttributeRowFromShapeIndex(i).getValue(length_col);
                });
        }
        for (auto& th : threads) th.join();
    }

    // --- Inicialização das variáveis ---
    int radiussize = radius.size();
    int radiusmask = 0;
    for (int i = 0; i < radiussize; i++) radiusmask |= (1 << i);

    // Acúmulo global final
    struct ChoiceAccum { double choice = 0, weighted_choice = 0, weighted_choice2 = 0; };
    // Cada thread terá sua própria cópia de global_choice_thread[thread_id][...]
    std::vector<std::vector<std::vector<std::array<ChoiceAccum, 2>>>> global_choice_thread(
        num_threads,
        std::vector<std::vector<std::array<ChoiceAccum, 2>>>(
            n, std::vector<std::array<ChoiceAccum, 2>>(radiussize)
        )
    );

    // --- Worker Principal (forward + backward) ---
    auto forward_and_backward_worker = [&](int thread_id, size_t start, size_t end) {
        // Local para acumuladores desta thread
        auto& local_choice = global_choice_thread[thread_id];
        // --- Loop das origens ---

        // --- 1. Audit trail LOCAL (apenas para esta origem) ---
        std::vector<std::vector<std::array<AuditTrailCell, 2>>> audittrail_local(
            n, std::vector<std::array<AuditTrailCell, 2>>(radiussize)
        );
        for (size_t origin = start; origin < end; ++origin) {
            
            for (size_t dest = 0; dest < n; ++dest)
                for (int k = 0; k < radiussize; ++k)
                    for (int d = 0; d < 2; ++d)
                        audittrail_local[dest][k][d].clearLine();

            // --- 2. Forward traversal (DepthmapX) ---
            std::vector<std::vector<SegmentData>> bins(tulip_bins);
            for (int k = 0; k < tulip_bins; ++k) bins[k].clear();
            std::vector<std::array<int, 2>> uncovered(n, { radiusmask, radiusmask });

            double rootseglength = map.getAttributeRowFromShapeIndex(origin).getValue(length_col);
            SegmentData segmentData(0, origin, SegmentRef(), 0, 0.5 * rootseglength, radiusmask);
            auto it = std::lower_bound(bins[0].begin(), bins[0].end(), segmentData);
            if (it == bins[0].end() || segmentData != *it) bins[0].insert(it, segmentData);

            int depthlevel = 0, opencount = 1;
            size_t currentbin = 0;
            while (opencount) {
                while (bins[currentbin].empty()) {
                    depthlevel++;
                    currentbin = (currentbin + 1) % tulip_bins;
                }
                SegmentData lineindex = bins[currentbin].back();
                bins[currentbin].pop_back();
                opencount--;

                int ref = lineindex.ref;
                int dir = (lineindex.dir == 1) ? 0 : 1;
                int coverage = lineindex.coverage & uncovered[ref][dir];
                if (coverage != 0) {
                    int rbin = 0, rbinbase;
                    if (lineindex.previous.ref != -1) {
                        uncovered[ref][dir] &= ~coverage;
                        while (((coverage >> rbin) & 0x1) == 0) rbin++;
                        rbinbase = rbin;
                        while (rbin < radiussize) {
                            if (((coverage >> rbin) & 0x1) == 1) {
                                audittrail_local[ref][rbin][dir].depth = depthlevel;
                                audittrail_local[ref][rbin][dir].previous = lineindex.previous;
                                audittrail_local[lineindex.previous.ref][rbin][(lineindex.previous.dir == 1) ? 0 : 1].leaf = false;
                            }
                            rbin++;
                        }
                    }
                    else {
                        rbinbase = 0;
                        uncovered[ref][0] &= ~coverage;
                        uncovered[ref][1] &= ~coverage;
                    }
                    Connector& line = map.getConnections()[ref];
                    float seglength;
                    int extradepth;
                    if (lineindex.dir != -1) {
                        for (auto& segconn : line.m_forward_segconns) {
                            rbin = rbinbase;
                            SegmentRef conn = segconn.first;
                            if ((uncovered[conn.ref][(conn.dir == 1 ? 0 : 1)] & coverage) != 0) {
                                extradepth = (routeweight_col != -1)
                                    ? (int)floor(segconn.second * tulip_bins * 0.5 * routeweights[conn.ref])
                                    : (int)floor(segconn.second * tulip_bins * 0.5);
                                seglength = lengths[conn.ref];
                                switch (m_radius_type) {
                                case Options::RADIUS_ANGULAR:
                                    while (rbin != radiussize && radius[rbin] != -1 && depthlevel + extradepth > (int)radius[rbin]) rbin++;
                                    break;
                                case Options::RADIUS_METRIC:
                                    while (rbin != radiussize && radius[rbin] != -1 && lineindex.metricdepth + seglength * 0.5 > radius[rbin]) rbin++;
                                    break;
                                case Options::RADIUS_STEPS:
                                    if (rbin != radiussize && radius[rbin] != -1 && lineindex.segdepth >= (int)radius[rbin]) rbin++;
                                    break;
                                }
                                if ((coverage >> rbin) != 0) {
                                    SegmentData sd(conn, SegmentRef(1, lineindex.ref), lineindex.segdepth + 1,
                                        lineindex.metricdepth + seglength, (coverage >> rbin) << rbin);
                                    size_t bin = (currentbin + tulip_bins + extradepth) % tulip_bins;
                                    depthmapX::insert_sorted(bins[bin], sd);
                                    opencount++;
                                }
                            }
                        }
                    }
                    if (lineindex.dir != 1) {
                        for (auto& segconn : line.m_back_segconns) {
                            rbin = rbinbase;
                            SegmentRef conn = segconn.first;
                            if ((uncovered[conn.ref][(conn.dir == 1 ? 0 : 1)] & coverage) != 0) {
                                extradepth = (routeweight_col != -1)
                                    ? (int)floor(segconn.second * tulip_bins * 0.5 * routeweights[conn.ref])
                                    : (int)floor(segconn.second * tulip_bins * 0.5);
                                seglength = lengths[conn.ref];
                                switch (m_radius_type) {
                                case Options::RADIUS_ANGULAR:
                                    while (rbin != radiussize && radius[rbin] != -1 && depthlevel + extradepth > (int)radius[rbin]) rbin++;
                                    break;
                                case Options::RADIUS_METRIC:
                                    while (rbin != radiussize && radius[rbin] != -1 && lineindex.metricdepth + seglength * 0.5 > radius[rbin]) rbin++;
                                    break;
                                case Options::RADIUS_STEPS:
                                    if (rbin != radiussize && radius[rbin] != -1 && lineindex.segdepth >= (int)radius[rbin]) rbin++;
                                    break;
                                }
                                if ((coverage >> rbin) != 0) {
                                    SegmentData sd(conn, SegmentRef(-1, lineindex.ref), lineindex.segdepth + 1,
                                        lineindex.metricdepth + seglength, (coverage >> rbin) << rbin);
                                    size_t bin = (currentbin + tulip_bins + extradepth) % tulip_bins;
                                    depthmapX::insert_sorted(bins[bin], sd);
                                    opencount++;
                                }
                            }
                        }
                    }
                }
            }

            // --- 3. Métricas integration/depth para essa origem (opcional, direto do audittrail_local) ---
            AttributeRow& row = attributes.getRow(AttributeKey(depthmapX::getMapAtIndex(map.getAllShapes(), origin)->first));
            for (int k = 0; k < radiussize; ++k) {
                double curs_node_count = 0.0, curs_total_depth = 0.0;
                double curs_total_weight = 0.0, curs_total_weighted_depth = 0.0;
                for (size_t dest = 0; dest < n; ++dest) {
                    for (int d = 0; d < 2; ++d) {
                        if (origin == dest) continue;
                        if (audittrail_local[dest][k][d].depth > 0) {
                            curs_node_count++;
                            curs_total_depth += audittrail_local[dest][k][d].depth;
                            curs_total_weight += weights[dest];
                            curs_total_weighted_depth += audittrail_local[dest][k][d].depth * weights[dest];
                        }
                    }
                }
                double total_depth_conv = curs_total_depth / ((tulip_bins - 1.0f) * 0.5f);
                double total_weighted_depth_conv = curs_total_weighted_depth / ((tulip_bins - 1.0f) * 0.5f);

                row.setValue(count_col[k], float(curs_node_count));
                if (curs_node_count > 1) {
                    row.setValue(td_col[k], total_depth_conv);
                    if (m_weighted_measure_col != -1) {
                        row.setValue(total_weight_col[k], float(curs_total_weight));
                        row.setValue(w_td_col[k], float(total_weighted_depth_conv));
                    }
                }
                else {
                    row.setValue(td_col[k], -1);
                    if (m_weighted_measure_col != -1) {
                        row.setValue(total_weight_col[k], -1.0f);
                        row.setValue(w_td_col[k], -1.0f);
                    }
                }
                if (total_depth_conv > 1e-9) {
                    row.setValue(integ_col[k], (float)(curs_node_count * curs_node_count / total_depth_conv));
                    if (m_weighted_measure_col != -1) {
                        row.setValue(w_integ_col[k], (float)(curs_total_weight * curs_total_weight / total_weighted_depth_conv));
                    }
                }
                else {
                    row.setValue(integ_col[k], -1);
                    if (m_weighted_measure_col != -1) {
                        row.setValue(w_integ_col[k], -1.0f);
                    }
                }
            }

            // --- 4. Backward: calcula o choice desta origem acumulando em local_choice ---
            // (O backward_worker é incorporado direto aqui, mas você pode chamar uma função se preferir)
            std::vector<std::vector<std::array<std::unordered_set<size_t>, 2>>> local_choicecovered(
                n, std::vector<std::array<std::unordered_set<size_t>, 2>>(radiussize)
            );

            for (int k = 0; k < radiussize; ++k) {
                for (size_t dest = 0; dest < n; ++dest) {
                    for (int d = 0; d < 2; ++d) {
                        if (audittrail_local[dest][k][d].leaf && (origin != dest)) {
                            int choicecount = 0;
                            double choiceweight = 0.0, choiceweight2 = 0.0;
                            SegmentRef here = SegmentRef(d == 0 ? 1 : -1, dest);
                            while (here.ref != static_cast<int>(origin)) {
                                int heredir = (here.dir == 1) ? 0 : 1;
                                if (here.ref < 0 || static_cast<size_t>(here.ref) >= n)
                                    break;
                                // USO DO UNORDERED_SET:
                                if (local_choicecovered[origin][k][heredir].count(here.ref) == 0) {
                                    choicecount++;
                                    choiceweight += weights[here.ref] * weights[origin];
                                    choiceweight2 += weights2[here.ref] * weights[origin];
                                    local_choicecovered[origin][k][heredir].insert(here.ref);
                                    if (m_weighted_measure_col != -1) {
                                        local_choice[here.ref][k][heredir].weighted_choice += (weights[here.ref] * weights[origin]) / 2.0;
                                        if (weighting_col2 != -1) {
                                            local_choice[here.ref][k][heredir].weighted_choice2 += (weights2[here.ref] * weights[origin]) / 2.0;
                                        }
                                    }
                                }
                                local_choice[here.ref][k][heredir].choice += choicecount;
                                local_choice[here.ref][k][heredir].weighted_choice += choiceweight;
                                local_choice[here.ref][k][heredir].weighted_choice2 += choiceweight2;
                                here = audittrail_local[here.ref][k][heredir].previous;
                            }
                            // No nó inicial (origin)
                            if (m_weighted_measure_col != -1 && here.ref == static_cast<int>(origin)) {
                                int origindir = (here.dir == 1) ? 0 : 1;
                                local_choice[here.ref][k][origindir].weighted_choice += choiceweight / 2.0;
                                if (weighting_col2 != -1) {
                                    local_choice[here.ref][k][origindir].weighted_choice2 += choiceweight2 / 2.0;
                                }
                            }
                        }
                    }
                }
            }

            // --- 5. Barra de progresso (opcional, pode remover para benchmark real) ---
            size_t value = ++global_processed_rows;
            if (comm) {
                if (qtimer(atime, 500)) {
                    if (comm->IsCancelled()) throw Communicator::CancelledException();
                    comm->CommPostMessage(Communicator::CURRENT_RECORD, value);
                }
            }
        } // (fim do loop de origin)
        };

    // --- Threading: cada thread processa seu chunk, armazena resultado no seu buffer local ---
    size_t chunk = (n + num_threads - 1) / num_threads;
    std::vector<std::thread> threads;
    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk;
        size_t end = std::min(n, (t + 1) * chunk);
        threads.emplace_back(forward_and_backward_worker, t, start, end);
    }
    for (auto& th : threads) th.join();

    // --- Redução Final: soma tudo em global_choice (sem lock, agora é rápido) ---
    std::vector<std::vector<std::array<ChoiceAccum, 2>>> global_choice(
        n, std::vector<std::array<ChoiceAccum, 2>>(radiussize)
    );
    for (int t = 0; t < num_threads; ++t)
        for (size_t i = 0; i < n; ++i)
            for (int r = 0; r < radiussize; ++r)
                for (int d = 0; d < 2; ++d) {
                    global_choice[i][r][d].choice += global_choice_thread[t][i][r][d].choice;
                    global_choice[i][r][d].weighted_choice += global_choice_thread[t][i][r][d].weighted_choice;
                    global_choice[i][r][d].weighted_choice2 += global_choice_thread[t][i][r][d].weighted_choice2;
                }


    // --- Atualiza atributos finais ---
    if (m_choice) {
        for (size_t cursor = 0; cursor < n; cursor++) {
            AttributeRow& row = attributes.getRow(AttributeKey(depthmapX::getMapAtIndex(map.getAllShapes(), cursor)->first));
            for (size_t r = 0; r < radiussize; r++) {
                double total_choice = global_choice[cursor][r][0].choice + global_choice[cursor][r][1].choice;
                double total_weighted_choice = global_choice[cursor][r][0].weighted_choice + global_choice[cursor][r][1].weighted_choice;
                double total_weighted_choice2 = global_choice[cursor][r][0].weighted_choice2 + global_choice[cursor][r][1].weighted_choice2;

                row.setValue(choice_col[r], float(total_choice));
                if (m_weighted_measure_col != -1) {
                    row.setValue(w_choice_col[r], float(total_weighted_choice));
                    if (weighting_col2 != -1) {
                        row.setValue(w_choice_col2[r], float(total_weighted_choice2));
                    }
                }
            }
        }
    }

    map.setDisplayedAttribute(-2);
    if (m_choice) map.setDisplayedAttribute(choice_col.back());
    else map.setDisplayedAttribute(td_col.back());
    return global_processed_rows > 0;
}
