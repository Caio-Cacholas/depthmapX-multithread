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

#include "salalib/axialmodules/axialintegration.h"

#include "genlib/pflipper.h"
#include <thread>
#include <vector>
#include <atomic>
#include <algorithm>
#include "genlib/stringutils.h"

bool AxialIntegration::run(Communicator* comm, ShapeGraph& map, bool simple_version) {
    // 1. Início: Setup de tempo e colunas
    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS, map.getShapeCount());
    }

    AttributeTable& attributes = map.getAttributeTable();

    // 2. Determinação dos raios
    bool radius_n = false;
    std::vector<int> radii;
    for (double radius : m_radius_set) {
        if (radius < 0) {
            radius_n = true;
        }
        else {
            radii.push_back(static_cast<int>(radius));
        }
    }
    if (radius_n) {
        radii.push_back(-1);
    }

    // 3. Pesos, se for análise ponderada
    std::vector<double> weights;
    std::string weighting_col_text;
    if (m_weighted_measure_col != -1) {
        weighting_col_text = attributes.getColumnName(m_weighted_measure_col);
        for (size_t i = 0; i < map.getShapeCount(); i++) {
            weights.push_back(map.getAttributeRowFromShapeIndex(i).getValue(m_weighted_measure_col));
        }
    }

    // 4. Adição de colunas necessárias
    for (int radius : radii) {
        std::string radius_text;
        if (radius != -1) radius_text = dXstring::formatString(radius, " R%d");
        if (m_choice) {
            attributes.insertOrResetColumn((std::string("Choice") + radius_text).c_str());
            attributes.insertOrResetColumn((std::string("Choice [Norm]") + radius_text).c_str());
            if (m_weighted_measure_col != -1) {
                attributes.insertOrResetColumn((std::string("Choice [") + weighting_col_text + " Wgt]" + radius_text).c_str());
                attributes.insertOrResetColumn((std::string("Choice [") + weighting_col_text + " Wgt][Norm]" + radius_text).c_str());
            }
        }
        if (!simple_version) {
            attributes.insertOrResetColumn((std::string("Entropy") + radius_text).c_str());
        }
        attributes.insertOrResetColumn((std::string("Integration [HH]") + radius_text).c_str());
        if (!simple_version) {
            attributes.insertOrResetColumn((std::string("Integration [P-value]") + radius_text).c_str());
            attributes.insertOrResetColumn((std::string("Integration [Tekl]") + radius_text).c_str());
            attributes.insertOrResetColumn((std::string("Intensity") + radius_text).c_str());
            attributes.insertOrResetColumn((std::string("Harmonic Mean Depth") + radius_text).c_str());
        }
        attributes.insertOrResetColumn((std::string("Mean Depth") + radius_text).c_str());
        attributes.insertOrResetColumn((std::string("Node Count") + radius_text).c_str());
        if (!simple_version) {
            attributes.insertOrResetColumn((std::string("Relativised Entropy") + radius_text).c_str());
        }
        if (m_weighted_measure_col != -1) {
            attributes.insertOrResetColumn((std::string("Mean Depth [") + weighting_col_text + " Wgt]" + radius_text).c_str());
            attributes.insertOrResetColumn((std::string("Total ") + weighting_col_text + radius_text).c_str());
        }
        if (m_fulloutput) {
            if (!simple_version) {
                attributes.insertOrResetColumn((std::string("RA [Penn]") + radius_text).c_str());
            }
            attributes.insertOrResetColumn((std::string("RA") + radius_text).c_str());
            if (!simple_version) {
                attributes.insertOrResetColumn((std::string("RRA") + radius_text).c_str());
            }
            attributes.insertOrResetColumn((std::string("Total Depth") + radius_text).c_str());
        }
    }
    if (m_local && !simple_version) {
        attributes.insertOrResetColumn("Control");
        attributes.insertOrResetColumn("Controllability");
    }

    // 5. Lookup dos índices das colunas
    std::vector<int> choice_col, n_choice_col, w_choice_col, nw_choice_col, entropy_col, integ_dv_col, integ_pv_col,
        integ_tk_col, intensity_col, depth_col, count_col, rel_entropy_col, penn_norm_col, w_depth_col,
        total_weight_col, ra_col, rra_col, td_col, harmonic_col;
    for (int radius : radii) {
        std::string radius_text;
        if (radius != -1) radius_text = std::string(" R") + dXstring::formatString(int(radius), "%d");
        if (m_choice) {
            choice_col.push_back(attributes.getColumnIndex((std::string("Choice") + radius_text).c_str()));
            n_choice_col.push_back(attributes.getColumnIndex((std::string("Choice [Norm]") + radius_text).c_str()));
            if (m_weighted_measure_col != -1) {
                w_choice_col.push_back(attributes.getColumnIndex((std::string("Choice [") + weighting_col_text + " Wgt]" + radius_text).c_str()));
                nw_choice_col.push_back(attributes.getColumnIndex((std::string("Choice [") + weighting_col_text + " Wgt][Norm]" + radius_text).c_str()));
            }
        }
        if (!simple_version) {
            entropy_col.push_back(attributes.getColumnIndex((std::string("Entropy") + radius_text).c_str()));
        }
        integ_dv_col.push_back(attributes.getColumnIndex((std::string("Integration [HH]") + radius_text).c_str()));
        if (!simple_version) {
            integ_pv_col.push_back(attributes.getColumnIndex((std::string("Integration [P-value]") + radius_text).c_str()));
            integ_tk_col.push_back(attributes.getColumnIndex((std::string("Integration [Tekl]") + radius_text).c_str()));
            intensity_col.push_back(attributes.getColumnIndex((std::string("Intensity") + radius_text).c_str()));
            harmonic_col.push_back(attributes.getColumnIndex((std::string("Harmonic Mean Depth") + radius_text).c_str()));
        }
        depth_col.push_back(attributes.getColumnIndex((std::string("Mean Depth") + radius_text).c_str()));
        count_col.push_back(attributes.getColumnIndex((std::string("Node Count") + radius_text).c_str()));
        if (!simple_version) {
            rel_entropy_col.push_back(attributes.getColumnIndex((std::string("Relativised Entropy") + radius_text).c_str()));
        }
        if (m_weighted_measure_col != -1) {
            w_depth_col.push_back(attributes.getColumnIndex((std::string("Mean Depth [") + weighting_col_text + " Wgt]" + radius_text).c_str()));
            total_weight_col.push_back(attributes.getColumnIndex((std::string("Total ") + weighting_col_text + radius_text).c_str()));
        }
        if (m_fulloutput) {
            ra_col.push_back(attributes.getColumnIndex((std::string("RA") + radius_text).c_str()));
            if (!simple_version) {
                penn_norm_col.push_back(attributes.getColumnIndex((std::string("RA [Penn]") + radius_text).c_str()));
                rra_col.push_back(attributes.getColumnIndex((std::string("RRA") + radius_text).c_str()));
            }
            td_col.push_back(attributes.getColumnIndex((std::string("Total Depth") + radius_text).c_str()));
        }
    }
    int control_col = -1, controllability_col = -1;
    if (m_local && !simple_version) {
        control_col = attributes.getColumnIndex("Control");
        controllability_col = attributes.getColumnIndex("Controllability");
    }

    // 6. PARTE MULTITHREAD
    size_t nShapes = map.getShapeCount();
    int num_threads = num_threadsGlobal;
    if (num_threads <= 0) num_threads = 1;
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (max_threads == 0) max_threads = 1;
    if (num_threads > static_cast<int>(max_threads)) num_threads = static_cast<int>(max_threads);
    size_t chunk = (nShapes + num_threads - 1) / num_threads;
    std::atomic<size_t> global_processed_rows(0);

    // --- Audittrail por thread (apenas se m_choice) ---
    std::vector<AnalysisInfo**> audittrail_per_thread(num_threads, nullptr);
    if (m_choice) {
        for (int t = 0; t < num_threads; ++t) {
            audittrail_per_thread[t] = new AnalysisInfo * [nShapes];
            for (size_t i = 0; i < nShapes; ++i)
                audittrail_per_thread[t][i] = new AnalysisInfo[radii.size()];
        }
    }

    auto process = [&](int thread_id, size_t start, size_t end) {
        time_t local_time = 0;
        bool* covered = new bool[nShapes];
        AnalysisInfo** audittrail = m_choice ? audittrail_per_thread[thread_id] : nullptr;
        for (size_t i = start; i < end && i < nShapes; ++i) {
            AttributeRow& row = map.getAttributeRowFromShapeIndex(i);
            for (size_t j = 0; j < nShapes; ++j) covered[j] = false;
            if (m_choice) for (size_t k = 0; k < nShapes; ++k) audittrail[k][0].previous.ref = -1;

            // --- m_local (Control/Controllability) ---
            if (m_local) {
                double control = 0.0;
                const std::vector<int>& connections = map.getConnections()[i].m_connections;
                std::vector<int> totalneighbourhood;
                for (int connection : connections) {
                    depthmapX::addIfNotExists(totalneighbourhood, connection);
                    int retro_size = 0;
                    auto& retconnectors = map.getConnections()[size_t(connection)].m_connections;
                    for (auto retconnector : retconnectors) {
                        retro_size++;
                        depthmapX::addIfNotExists(totalneighbourhood, retconnector);
                    }
                    control += 1.0 / double(retro_size);
                }
                if (!simple_version) {
                    if (connections.size() > 0) {
                        row.setValue(control_col, float(control));
                        row.setValue(controllability_col,
                            float(double(connections.size()) / double(totalneighbourhood.size() - 1)));
                    }
                    else {
                        row.setValue(control_col, -1);
                        row.setValue(controllability_col, -1);
                    }
                }
            }

            std::vector<int> depthcounts; depthcounts.push_back(0);
            pflipper<std::vector<std::pair<int, int>>> foundlist;
            foundlist.a().push_back(std::pair<int, int>(i, -1));
            covered[i] = true;
            int total_depth = 0, depth = 1, node_count = 1, pos = -1, previous = -1;
            double weight = 0.0, rootweight = 0.0, total_weight = 0.0, w_total_depth = 0.0;
            if (m_weighted_measure_col != -1) {
                rootweight = weights[i];
                total_weight += rootweight;
            }
            int index = -1;
            int r = 0;
            for (int radius : radii) {
                while (foundlist.a().size()) {
                    if (!m_choice) {
                        index = foundlist.a().back().first;
                    }
                    else {
                        pos = pafrand() % foundlist.a().size();
                        index = foundlist.a().at(pos).first;
                        previous = foundlist.a().at(pos).second;
                        audittrail[index][0].previous.ref = previous;
                    }
                    Connector& line = map.getConnections()[index];
                    for (size_t k = 0; k < line.m_connections.size(); ++k) {
                        if (!covered[line.m_connections[k]]) {
                            covered[line.m_connections[k]] = true;
                            foundlist.b().push_back(std::pair<int, int>(line.m_connections[k], index));
                            if (m_weighted_measure_col != -1) {
                                weight = weights[line.m_connections[k]];
                                total_weight += weight;
                                w_total_depth += depth * weight;
                            }
                            if (m_choice && previous != -1) {
                                size_t here = index;
                                while (here != i) {
                                    if (audittrail[here][0].previous.ref == -1) break;
                                    audittrail[here][r].choice += 1;
                                    audittrail[here][r].weighted_choice += weight * rootweight;
                                    here = audittrail[here][0].previous.ref;
                                }
                                if (m_weighted_measure_col != -1) {
                                    audittrail[i][r].weighted_choice += (weight * rootweight) * 0.5;
                                    audittrail[line.m_connections[k]][r].weighted_choice += (weight * rootweight) * 0.5;
                                }
                            }
                            total_depth += depth;
                            node_count++;
                            depthcounts.back() += 1;
                        }
                    }
                    if (!m_choice)
                        foundlist.a().pop_back();
                    else
                        foundlist.a().erase(foundlist.a().begin() + pos);
                    if (!foundlist.a().size()) {
                        foundlist.flip();
                        depth++;
                        depthcounts.push_back(0);
                        if (radius != -1 && depth > radius) break;
                    }
                }
                row.setValue(count_col[r], float(node_count));
                if (m_weighted_measure_col != -1) row.setValue(total_weight_col[r], float(total_weight));
                if (node_count > 1) {
                    double mean_depth = double(total_depth) / double(node_count - 1);
                    row.setValue(depth_col[r], float(mean_depth));
                    if (m_weighted_measure_col != -1) row.setValue(w_depth_col[r], float(w_total_depth / total_weight));
                    if (node_count > 2 && mean_depth > 1.0) {
                        double ra = 2.0 * (mean_depth - 1.0) / double(node_count - 2);
                        double rra_d = ra / dvalue(node_count);
                        double rra_p = ra / dvalue(node_count);
                        double integ_tk = teklinteg(node_count, total_depth);
                        row.setValue(integ_dv_col[r], float(1.0 / rra_d));
                        if (!simple_version) {
                            row.setValue(integ_pv_col[r], float(1.0 / rra_p));
                            if (total_depth - node_count + 1 > 1)
                                row.setValue(integ_tk_col[r], float(integ_tk));
                            else
                                row.setValue(integ_tk_col[r], -1.0f);
                        }
                        if (m_fulloutput) {
                            row.setValue(ra_col[r], float(ra));
                            if (!simple_version) row.setValue(rra_col[r], float(rra_d));
                            row.setValue(td_col[r], float(total_depth));
                            if (!simple_version) {
                                double dmin = node_count - 1;
                                double dmax = palmtree(node_count, depth - 1);
                                if (dmax != dmin)
                                    row.setValue(penn_norm_col[r], float((dmax - total_depth) / (dmax - dmin)));
                            }
                        }
                    }
                    else {
                        row.setValue(integ_dv_col[r], -1.0f);
                        if (!simple_version) {
                            row.setValue(integ_pv_col[r], -1.0f);
                            row.setValue(integ_tk_col[r], -1.0f);
                        }
                        if (m_fulloutput) {
                            row.setValue(ra_col[r], -1.0f);
                            if (!simple_version) row.setValue(rra_col[r], -1.0f);
                            row.setValue(td_col[r], -1.0f);
                            if (!simple_version) row.setValue(penn_norm_col[r], -1.0f);
                        }
                    }
                    if (!simple_version) {
                        double entropy = 0.0, intensity = 0.0, rel_entropy = 0.0, factorial = 1.0, harmonic = 0.0;
                        for (size_t k = 0; k < depthcounts.size(); k++) {
                            if (depthcounts[k] != 0) {
                                double prob = double(depthcounts[k]) / double(node_count);
                                entropy -= prob * log2(prob);
                                factorial *= double(k + 1);
                                double q = (pow(mean_depth, double(k)) / double(factorial)) * exp(-mean_depth);
                                rel_entropy += (double)prob * log2(prob / q);
                                harmonic += 1.0 / double(depthcounts[k]);
                            }
                        }
                        harmonic = double(depthcounts.size()) / harmonic;
                        if (total_depth > node_count)
                            intensity = node_count * entropy / (total_depth - node_count);
                        else
                            intensity = -1;
                        row.setValue(entropy_col[r], float(entropy));
                        row.setValue(rel_entropy_col[r], float(rel_entropy));
                        row.setValue(intensity_col[r], float(intensity));
                        row.setValue(harmonic_col[r], float(harmonic));
                    }
                }
                else {
                    row.setValue(depth_col[r], -1.0f);
                    row.setValue(integ_dv_col[r], -1.0f);
                    if (!simple_version) {
                        row.setValue(integ_pv_col[r], -1.0f);
                        row.setValue(integ_tk_col[r], -1.0f);
                        row.setValue(entropy_col[r], -1.0f);
                        row.setValue(rel_entropy_col[r], -1.0f);
                        row.setValue(harmonic_col[r], -1.0f);
                    }
                }
                ++r;
            }

            size_t processed = ++global_processed_rows;
            if (comm) {
                if (qtimer(local_time, 500)) {
                    if (comm->IsCancelled()) {
                        delete[] covered;
                        throw Communicator::CancelledException();
                    }
                    comm->CommPostMessage(Communicator::CURRENT_RECORD, processed);
                }
            }
        }

        delete[] covered;
        }; // end process lambda

    // 7. Criação e join das threads
    std::vector<std::thread> threads;
    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk;
        size_t end = std::min(nShapes, (t + 1) * chunk);
        threads.emplace_back(process, t, start, end);
    }
    for (auto& th : threads) th.join();

    // 8. Pós-processamento: Merge audittrail (caso m_choice)
    if (m_choice) {
        AnalysisInfo** audittrail_global = new AnalysisInfo * [nShapes];
        for (size_t i = 0; i < nShapes; ++i) {
            audittrail_global[i] = new AnalysisInfo[radii.size()];
            for (size_t r = 0; r < radii.size(); ++r)
                audittrail_global[i][r] = AnalysisInfo();
        }
        for (int t = 0; t < num_threads; ++t) {
            for (size_t i = 0; i < nShapes; ++i)
                for (size_t r = 0; r < radii.size(); ++r) {
                    audittrail_global[i][r].choice += audittrail_per_thread[t][i][r].choice;
                    audittrail_global[i][r].weighted_choice += audittrail_per_thread[t][i][r].weighted_choice;
                }
        }

        // Preenchimento dos valores nas colunas de choice
        size_t i = -1;
        for (auto& iter : attributes) {
            i++;
            AttributeRow& row = iter.getRow();
            double total_choice = 0.0, w_total_choice = 0.0;
            for (size_t r = 0; r < radii.size(); r++) {
                total_choice += audittrail_global[i][r].choice;
                w_total_choice += audittrail_global[i][r].weighted_choice;
                double node_count = row.getValue(count_col[r]);
                double total_weight = 0;
                if (m_weighted_measure_col != -1) {
                    total_weight = row.getValue(total_weight_col[r]);
                }
                if (node_count > 2) {
                    row.setValue(choice_col[r], float(total_choice));
                    row.setValue(n_choice_col[r], float(2.0 * total_choice / ((node_count - 1) * (node_count - 2))));
                    if (m_weighted_measure_col != -1) {
                        row.setValue(w_choice_col[r], float(w_total_choice));
                        row.setValue(nw_choice_col[r], float(2.0 * w_total_choice / (total_weight * total_weight)));
                    }
                }
                else {
                    row.setValue(choice_col[r], -1);
                    row.setValue(n_choice_col[r], -1);
                    if (m_weighted_measure_col != -1) {
                        row.setValue(w_choice_col[r], -1);
                        row.setValue(nw_choice_col[r], -1);
                    }
                }
            }
        }
        for (size_t i = 0; i < nShapes; ++i) delete[] audittrail_global[i];
        delete[] audittrail_global;
        for (int t = 0; t < num_threads; ++t) {
            for (size_t i = 0; i < nShapes; ++i) delete[] audittrail_per_thread[t][i];
            delete[] audittrail_per_thread[t];
        }
    }

    // 9. Display
    map.setDisplayedAttribute(-1);
    map.setDisplayedAttribute(integ_dv_col.back());
    return true;
}
