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

#include "salalib/segmmodules/segmangular.h"
#include "salalib/options.h"
#include "genlib/stringutils.h"
#include <thread>
#include <atomic>
#include <vector>
#include <cmath>
#include <algorithm>

// Define struct para armazenar os resultados intermediários por shape/radius
struct AngularResult {
    std::vector<float> count;
    std::vector<float> depth;
    std::vector<float> total;
    AngularResult(size_t nradii) : count(nradii, 0), depth(nradii, -1), total(nradii, -1) {}
};

bool SegmentAngular::run(Communicator *comm, ShapeGraph &map, bool) {

    if (map.getMapType() != ShapeMap::SEGMENTMAP) {
        return false;
    }

    AttributeTable &attributes = map.getAttributeTable();

    // >>> CONTADOR ATÔMICO GLOBAL <<<
    std::atomic<size_t> global_processed_rows(0);

    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS, map.getConnections().size());
    }

    // --- Calcula os bins de raio ---
    bool radius_n = false;
    std::vector<double> radii;
    for (double radius : m_radius_set) {
        if (radius < 0) {
            radius_n = true;
        } else {
            radii.push_back(radius);
        }
    }
    if (radius_n) {
        radii.push_back(-1.0);
    }

    // --- Insere colunas na tabela de atributos ---
    std::vector<int> depth_col, count_col, total_col;
    for (double radius : radii) {
        std::string radius_text = makeRadiusText(Options::RADIUS_ANGULAR, radius);
        std::string depth_col_text = std::string("Angular Mean Depth") + radius_text;
        attributes.insertOrResetColumn(depth_col_text.c_str());
        std::string count_col_text = std::string("Angular Node Count") + radius_text;
        attributes.insertOrResetColumn(count_col_text.c_str());
        std::string total_col_text = std::string("Angular Total Depth") + radius_text;
        attributes.insertOrResetColumn(total_col_text.c_str());
    }
    for (double radius : radii) {
        std::string radius_text = makeRadiusText(Options::RADIUS_ANGULAR, radius);
        depth_col.push_back(attributes.getColumnIndex((std::string("Angular Mean Depth") + radius_text).c_str()));
        count_col.push_back(attributes.getColumnIndex((std::string("Angular Node Count") + radius_text).c_str()));
        total_col.push_back(attributes.getColumnIndex((std::string("Angular Total Depth") + radius_text).c_str()));
    }

    // --- Multithreading ---
    size_t nShapes = map.getShapeCount();
    size_t nRadii = radii.size();

    // Vetor de resultados intermediários (um por shape)
    std::vector<AngularResult> results(nShapes, AngularResult(nRadii));

    int num_threads = num_threadsGlobal;
    if (num_threads <= 0) num_threads = 1;
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (max_threads == 0) max_threads = 1;
    if (num_threads > static_cast<int>(max_threads)) num_threads = static_cast<int>(max_threads);

    auto process = [&](size_t start, size_t end) {
        std::vector<bool> covered(nShapes);
        time_t local_time = 0;  // Para qtimer local por thread
        for (size_t i = start; i < end; ++i) {
            for (size_t j = 0; j < nShapes; j++) covered[j] = false;

            std::vector<std::pair<float, SegmentData>> anglebins;
            anglebins.push_back(std::make_pair(0.0f, SegmentData(0, i, SegmentRef(), 0, 0.0, 0)));

            std::vector<double> total_depth(nRadii, 0.0);
            std::vector<int> node_count(nRadii, 0);

            // --- Algoritmo principal ---
            while (anglebins.size()) {
                auto iter = anglebins.begin();
                SegmentData lineindex = iter->second;
                if (!covered[lineindex.ref]) {
                    covered[lineindex.ref] = true;
                    double depth_to_line = iter->first;
                    total_depth[lineindex.coverage] += depth_to_line;
                    node_count[lineindex.coverage] += 1;
                    anglebins.erase(iter);
                    Connector &line = map.getConnections()[lineindex.ref];
                    if (lineindex.dir != -1) {
                        for (auto &segconn : line.m_forward_segconns) {
                            if (!covered[segconn.first.ref]) {
                                double angle = depth_to_line + segconn.second;
                                size_t rbin = lineindex.coverage;
                                while (rbin != nRadii && radii[rbin] != -1 && angle > radii[rbin]) {
                                    rbin++;
                                }
                                if (rbin != nRadii) {
                                    depthmapX::insert_sorted(
                                        anglebins, std::make_pair(float(angle),
                                            SegmentData(segconn.first, SegmentRef(), 0, 0.0, rbin)));
                                }
                            }
                        }
                    }
                    if (lineindex.dir != 1) {
                        for (auto &segconn : line.m_back_segconns) {
                            if (!covered[segconn.first.ref]) {
                                double angle = depth_to_line + segconn.second;
                                size_t rbin = lineindex.coverage;
                                while (rbin != nRadii && radii[rbin] != -1 && angle > radii[rbin]) {
                                    rbin++;
                                }
                                if (rbin != nRadii) {
                                    depthmapX::insert_sorted(
                                        anglebins, std::make_pair(float(angle),
                                            SegmentData(segconn.first, SegmentRef(), 0, 0.0, rbin)));
                                }
                            }
                        }
                    }
                } else {
                    anglebins.erase(iter);
                }
            }
            // --- Finaliza atributos acumulados para o shape i ---
            int curs_node_count = 0;
            double curs_total_depth = 0.0;
            for (size_t r = 0; r < nRadii; r++) {
                curs_node_count += node_count[r];
                curs_total_depth += total_depth[r];
                results[i].count[r] = float(curs_node_count);
                if (curs_node_count > 1) {
                    results[i].depth[r] = float(curs_total_depth / double(curs_node_count - 1));
                    results[i].total[r] = float(curs_total_depth);
                } else {
                    results[i].depth[r] = -1;
                    results[i].total[r] = -1;
                }
            }
            // >>> ATUALIZAÇÃO DA BARRA DE PROGRESSO (dentro das threads) <<<
            size_t processed = ++global_processed_rows;  // incrementa e pega valor atual
            if (comm) {
                if (qtimer(local_time, 500)) { // só posta a cada 500ms por thread
                    if (comm->IsCancelled()) throw Communicator::CancelledException();
                    comm->CommPostMessage(Communicator::CURRENT_RECORD, processed);
                }
            }
        }
    };

    // Lança threads
    std::vector<std::thread> threads;
    size_t chunk = nShapes / num_threads;
    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk;
        size_t end = (t == num_threads - 1) ? nShapes : start + chunk;
        threads.emplace_back(process, start, end);
    }
    for (auto &th : threads) th.join();

    // --- Atualiza os atributos na tabela ---
    auto it = attributes.begin();
    for (size_t i = 0; i < nShapes; ++i, ++it) {
        AttributeRow& row = it->getRow();
        for (size_t r = 0; r < nRadii; ++r) {
            row.setValue(count_col[r], results[i].count[r]);
            row.setValue(depth_col[r], results[i].depth[r]);
            row.setValue(total_col[r], results[i].total[r]);
        }
        if (comm) {
            if (qtimer(atime, 500)) {
                if (comm->IsCancelled()) throw Communicator::CancelledException();
                comm->CommPostMessage(Communicator::CURRENT_RECORD, i);
            }
        }
    }

    map.setDisplayedAttribute(-2);
    map.setDisplayedAttribute(depth_col.back());

    return true;
}
