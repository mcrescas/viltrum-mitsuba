/*
 *  OpenVDB Data Source for Mitsuba Render
 *  Copyright (C) 2014 Bo Zhou<Bo.Schwarzstein@gmail.com>

 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <mitsuba/core/properties.h>
#include <mitsuba/render/volume.h>

MTS_NAMESPACE_BEGIN

class MTS_EXPORT_RENDER VdbDataSource : public VolumeDataSource
{
public :

    VdbDataSource(const Properties &props);

    VdbDataSource(Stream *stream, InstanceManager *manager);

    virtual ~VdbDataSource();

	virtual bool supportsFloatLookups() const;

	virtual Float lookupFloat(const Point &_p) const;

	virtual bool supportsVectorLookups() const;

    virtual bool supportsSpectrumLookups() const;

    virtual void serialize(Stream *stream, InstanceManager *manager) const;

	virtual Float getStepSize() const;

	virtual Float getMaximumFloatValue() const;

    virtual Float getAvgFloatValue() const;

    MTS_DECLARE_CLASS()

private :

    bool open();
    void getMinMaxVDB();
    void getAvgVDB();

    Float minVal;
    Float maxVal;
    Float avgVal;

    std::string m_filename;
    std::string m_fieldname;
    Float m_customStepSize;

    Transform m_worldToGrid;
    Transform m_worldToVolume;
    Transform m_volumeToWorld;
} ;

MTS_IMPLEMENT_CLASS_S(VdbDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(VdbDataSource, "OpenVDB data source");
MTS_NAMESPACE_END