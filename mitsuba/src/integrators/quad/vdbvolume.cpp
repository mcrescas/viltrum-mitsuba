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


#include <openvdb/openvdb.h>
#include <openvdb/math/BBox.h>
#include <openvdb/math/Ray.h>
#include <openvdb/tools/Interpolation.h>

#include "vdbvolume.h"

struct VdbGridSet
{
    std::map<std::string, openvdb::GridBase::Ptr > m_grids;
};

struct VdbGridPool
{
    std::map<std::string, VdbGridSet > m_gridSets;
};

static VdbGridPool s_vdbGridPool ;

MTS_NAMESPACE_BEGIN

VdbDataSource::VdbDataSource(const Properties &props)
    :
    VolumeDataSource(props),
    m_customStepSize(0.0f) {
    openvdb::initialize();

    //
    m_filename = props.getString("filename");
    m_fieldname = props.getString("fieldname");
    if( props.hasProperty("customStepSize") ) {
        m_customStepSize = props.getFloat("customStepSize");
    }

    std::map<std::string, VdbGridSet>::iterator gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);
    if (gridSetItr == s_vdbGridPool.m_gridSets.end()) {
        if (open()) {
            gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);

            std::map<std::string, openvdb::GridBase::Ptr >::iterator gridItr = gridSetItr->second.m_grids.find(m_fieldname);
            if (gridItr == gridSetItr->second.m_grids.end()) {
                Log(EError, "Can't get the specific field [%s] from [%s] to read.", m_fieldname.c_str(), m_filename.c_str());
            }
        }
        else {
            Log(EError, "Can't open the file [%s].", m_filename.c_str());
        }
    }
    else if (gridSetItr->second.m_grids.find(m_fieldname) == gridSetItr->second.m_grids.end()) {
        Log(EError, "Opened the file [%s] but can't get the field [%s].", m_filename.c_str(), m_fieldname.c_str());
    }

    /* Load transform */

    mitsuba::AABB m_dataAABB;
    Vector m_res;
    gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);
    if (gridSetItr != s_vdbGridPool.m_gridSets.end()) {
        std::map<std::string, openvdb::GridBase::Ptr >::const_iterator gridItr = gridSetItr->second.m_grids.find(m_fieldname);
        if (gridItr != gridSetItr->second.m_grids.end()) {
            openvdb::math::Coord voxelRes = gridItr->second->evalActiveVoxelDim();
            m_res[0] = voxelRes[0];
            m_res[1] = voxelRes[1];
            m_res[2] = voxelRes[2];

            openvdb::math::CoordBBox boundingBox = gridItr->second->evalActiveVoxelBoundingBox();
            Point minBBOX(boundingBox.min()[0], boundingBox.min()[1], boundingBox.min()[2]);
            Point maxBBOX(boundingBox.max()[0], boundingBox.max()[1], boundingBox.max()[2]);
            m_dataAABB = AABB(minBBOX, maxBBOX);
            Log(EInfo, "AABB data %s", m_dataAABB.toString().c_str());
        }
    }

    m_volumeToWorld = props.getTransform("toWorld", Transform());

    Vector extents(m_dataAABB.getExtents());
    m_worldToVolume = m_volumeToWorld.inverse();
    m_worldToGrid = Transform::scale(Vector(
            (m_res[0] - 1) / extents[0],
            (m_res[1] - 1) / extents[1],
            (m_res[2] - 1) / extents[2])
        ) * Transform::translate(-Vector(m_dataAABB.min)) * m_worldToVolume;

    // Precompute min/max values
    getMinMaxVDB();
    getAvgVDB();
}

VdbDataSource::VdbDataSource(Stream *stream, InstanceManager *manager)
    :
    VolumeDataSource(stream, manager) {
    m_filename = stream->readString();
    m_fieldname = stream->readString();
    m_customStepSize = stream->readFloat();
}

VdbDataSource::~VdbDataSource() {
}

bool VdbDataSource::supportsFloatLookups() const {
    return true;
}

Float VdbDataSource::lookupFloat(const Point &_p) const {
    //const Point p = m_worldToGrid.transformAffine(_p);
    const Point p = _p;
    //const Point p = m_worldToGrid(_p);
    Float data = 0.0f;

    std::map<std::string, VdbGridSet >::const_iterator gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);
    if (gridSetItr != s_vdbGridPool.m_gridSets.end()) {
        std::map<std::string, openvdb::GridBase::Ptr >::const_iterator gridItr = gridSetItr->second.m_grids.find(m_fieldname);
        if (gridItr != gridSetItr->second.m_grids.end()) {
            const openvdb::FloatGrid::Ptr &floatGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(gridItr->second);
            openvdb::tools::GridSampler<openvdb::FloatTree, openvdb::tools::BoxSampler> interpolator(floatGrid->constTree() , floatGrid->transform());
            data = interpolator.wsSample(openvdb::math::Vec3d(p.x, p.y, p.z));
        }
    }

    return data;
}

bool VdbDataSource::supportsVectorLookups() const {
    return false;
}

bool VdbDataSource::supportsSpectrumLookups() const {
    return false;
}

void VdbDataSource::serialize(Stream *stream, InstanceManager *manager) const {
    VolumeDataSource::serialize(stream, manager);

    stream->writeString(m_filename);
    stream->writeString(m_fieldname);
    stream->writeFloat(m_customStepSize);
}

Float VdbDataSource::getStepSize() const {
    Float stepSize = 1.0f;
    if ( m_customStepSize > 0.0f ) {
        stepSize = m_customStepSize;
    }
    std::map<std::string, VdbGridSet >::const_iterator gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);
    if (gridSetItr != s_vdbGridPool.m_gridSets.end()) {
        std::map<std::string, openvdb::GridBase::Ptr >::const_iterator gridItr = gridSetItr->second.m_grids.find(m_fieldname);
        if (gridItr != gridSetItr->second.m_grids.end()) {
            const openvdb::math::Vec3d &voxelSize = gridItr->second->constTransform().voxelSize();
            stepSize = std::min(std::min(voxelSize.x(), voxelSize.y()) , voxelSize.z()) * 0.5;
        }
    }
    return stepSize;
}

Float VdbDataSource::getMaximumFloatValue() const {
    return maxVal;
}

void VdbDataSource::getMinMaxVDB() {
    openvdb::FloatGrid::ValueType minV, maxV;
    std::map<std::string, VdbGridSet >::const_iterator gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);
    if (gridSetItr != s_vdbGridPool.m_gridSets.end()) {
        std::map<std::string, openvdb::GridBase::Ptr >::const_iterator gridItr = gridSetItr->second.m_grids.find(m_fieldname);
        if (gridItr != gridSetItr->second.m_grids.end()) {
            const openvdb::FloatGrid::Ptr &floatGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(gridItr->second);
            floatGrid->evalMinMax(minV, maxV);
            this->minVal = minV;
            this->maxVal = maxV;
        }
    }
}

Float VdbDataSource::getAvgFloatValue() const {
    return avgVal;
}

void VdbDataSource::getAvgVDB() {
    std::map<std::string, VdbGridSet >::const_iterator gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);
    if (gridSetItr != s_vdbGridPool.m_gridSets.end()) {
        std::map<std::string, openvdb::GridBase::Ptr >::const_iterator gridItr = gridSetItr->second.m_grids.find(m_fieldname);
        if (gridItr != gridSetItr->second.m_grids.end()) {
            const openvdb::FloatGrid::Ptr &floatGrid = openvdb::gridPtrCast<openvdb::FloatGrid>(gridItr->second);

            auto& tree = floatGrid->constTree();
            avgVal = 0.0;
            double count = 0;
            if (auto iter = tree.cbeginValueOn()) {
                avgVal = *iter;
                count++;
                for (++iter; iter; ++iter) {
                    const auto& val = *iter;
                    avgVal += val;
                    count++;
                }
            }
            avgVal /= count;
        }
    }
}

bool VdbDataSource::open() {
    bool success = false ;

    try {
        boost::scoped_ptr< openvdb::io::File > file(new openvdb::io::File(m_filename));
        if (file->open()) {
            s_vdbGridPool.m_gridSets.insert(std::make_pair(m_filename, VdbGridSet()));
            std::map<std::string, VdbGridSet>::iterator gridSetItr = s_vdbGridPool.m_gridSets.find(m_filename);

            //
            openvdb::math::Vec3d bbMin = openvdb::math::Vec3d(  std::numeric_limits< double >::max());
            openvdb::math::Vec3d bbMax = openvdb::math::Vec3d(- std::numeric_limits< double >::max());
            openvdb::GridPtrVecPtr grids = file->readAllGridMetadata() ;
            for (openvdb::GridPtrVec::const_iterator itr = grids->begin(); itr != grids->end(); ++ itr) {
                const openvdb::Vec3i & bbMinI = (*itr)->metaValue<openvdb::Vec3i>(openvdb::GridBase::META_FILE_BBOX_MIN);
                const openvdb::Vec3i & bbMaxI = (*itr)->metaValue<openvdb::Vec3i>(openvdb::GridBase::META_FILE_BBOX_MAX);
                bbMin = openvdb::math::minComponent(bbMin , (*itr)->indexToWorld(bbMinI));
                bbMax = openvdb::math::maxComponent(bbMax , (*itr)->indexToWorld(bbMaxI));
                const std::string &fieldname = (*itr)->getName();
                gridSetItr->second.m_grids.insert(std::make_pair(fieldname, file->readGrid((*itr)->getName())));
            }
            m_aabb = AABB(Point(bbMin.x(), bbMin.y(), bbMin.z()), Point(bbMax.x(), bbMax.y(), bbMax.z()));

            Log(EInfo, "AABB %s", m_aabb.toString().c_str());

            success = true;
        }
    }
    catch(const std::exception &e) {
        std::cerr << e.what() << std::endl;
    }
    catch(...) {
    }

    return success ;
}

MTS_NAMESPACE_END
