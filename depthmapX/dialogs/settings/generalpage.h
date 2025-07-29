// Copyright (C) 2017 Petros Koutsolampros

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

#pragma once

#include <QWidget>
#include <QSpinBox>
#include "thread_globals.h"
//#include "salalib/segmmodules/segmtopological.h"
//#include "salalib/segmmodules/segmtulip.h"
#include "settingspage.h"
#include <iostream>

class GeneralPage : public SettingsPage
{
private:
    bool m_simpleVersion = false;
    QSpinBox *m_spinThreads = nullptr;
    void readSettings(Settings &settings) {
        m_simpleVersion = settings.readSetting(SettingTag::simpleVersion, true).toBool();
        // Lê o valor de numThreads, default = 1
        int threads = settings.readSetting(SettingTag::numThreads, 1).toInt();
        if (m_spinThreads) m_spinThreads->setValue(threads);
        num_threadsGlobal = settings.readSetting(SettingTag::numThreads, 1).toInt();

    }
public:
    GeneralPage(Settings &settings, QWidget *parent = 0);
    virtual void writeSettings(Settings &settings) override {
        settings.writeSetting(SettingTag::simpleVersion, m_simpleVersion);
        // Salva o valor do spinbox (se já estiver inicializado)
        if (m_spinThreads) settings.writeSetting(SettingTag::numThreads, m_spinThreads->value());

    }
};
