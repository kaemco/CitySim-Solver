#include "voxel.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "RENIndexedFaceSet.h"
#include "RENIndexedFaceSetBuilder.h"
#include "building.h"
#include "climate.h"
#include "district.h"
#include "scene.h"
#include "surface.h"
#include "zone.h"

using namespace std;

namespace {

const uint8_t AIR_VALUE = static_cast<uint8_t>(VoxelMaterial::AIR);
const uint8_t BUILDING_VALUE = static_cast<uint8_t>(VoxelMaterial::BUILDING);
const uint8_t GROUND_VALUE = static_cast<uint8_t>(VoxelMaterial::GROUND);
const uint8_t TREE_VALUE = static_cast<uint8_t>(VoxelMaterial::TREE);
const uint8_t OUTSIDE_AIR_VALUE = 4;

struct Vec3 {
    Vec3():x(0.),y(0.),z(0.) {}
    Vec3(double xValue,double yValue,double zValue)
        :x(xValue),y(yValue),z(zValue) {}
    double x, y, z;
};

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return Vec3{a.x-b.x,a.y-b.y,a.z-b.z};
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y*b.z-a.z*b.y,
                a.z*b.x-a.x*b.z,
                a.x*b.y-a.y*b.x};
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

struct Triangle {
    Vec3 p[3];
    Vec3 minimum;
    Vec3 maximum;
};

Triangle makeTriangle(const Vec3& a, const Vec3& b, const Vec3& c) {
    Triangle triangle;
    triangle.p[0]=a; triangle.p[1]=b; triangle.p[2]=c;
    triangle.minimum = Vec3{min(a.x,min(b.x,c.x)),
                            min(a.y,min(b.y,c.y)),
                            min(a.z,min(b.z,c.z))};
    triangle.maximum = Vec3{max(a.x,max(b.x,c.x)),
                            max(a.y,max(b.y,c.y)),
                            max(a.z,max(b.z,c.z))};
    return triangle;
}

bool isDegenerate(const Triangle& triangle) {
    const Vec3 e0=triangle.p[1]-triangle.p[0];
    const Vec3 e1=triangle.p[2]-triangle.p[0];
    const Vec3 n=cross(e0,e1);
    return dot(n,n) <= 1e-20;
}

bool separatedOnAxis(const Vec3& axis, const Vec3 vertices[3], const Vec3& half) {
    const double lengthSquared=dot(axis,axis);
    if (lengthSquared <= 1e-24) return false;
    const double p0=dot(vertices[0],axis);
    const double p1=dot(vertices[1],axis);
    const double p2=dot(vertices[2],axis);
    const double minimum=min(p0,min(p1,p2));
    const double maximum=max(p0,max(p1,p2));
    const double radius=half.x*fabs(axis.x)
                       +half.y*fabs(axis.y)
                       +half.z*fabs(axis.z);
    const double tolerance=1e-11*(1.+radius);
    return minimum > radius+tolerance || maximum < -radius-tolerance;
}

/**
 * Generic 13-axis separating-axis triangle/AABB test.  It is slightly more
 * arithmetic than a specialized implementation but compact, branch-light and
 * robust for the anisotropic voxels used here.
 */
bool triangleIntersectsBox(const Triangle& triangle, const Vec3& centre, const Vec3& half) {
    Vec3 vertices[3] = {
        triangle.p[0]-centre,
        triangle.p[1]-centre,
        triangle.p[2]-centre
    };
    const Vec3 boxAxes[3] = {
        Vec3{1.,0.,0.}, Vec3{0.,1.,0.}, Vec3{0.,0.,1.}
    };
    for (unsigned int axis=0;axis<3;++axis)
        if (separatedOnAxis(boxAxes[axis],vertices,half)) return false;

    const Vec3 edges[3] = {
        vertices[1]-vertices[0],
        vertices[2]-vertices[1],
        vertices[0]-vertices[2]
    };
    const Vec3 normal=cross(edges[0],edges[1]);
    if (separatedOnAxis(normal,vertices,half)) return false;

    for (unsigned int edge=0;edge<3;++edge)
        for (unsigned int axis=0;axis<3;++axis)
            if (separatedOnAxis(cross(edges[edge],boxAxes[axis]),vertices,half))
                return false;
    return true;
}

struct IndexRange {
    size_t minimum=0;
    size_t maximum=0;
    bool valid=false;
};

IndexRange coordinateRange(double minimum, double maximum,
                           double origin, double delta, size_t count) {
    IndexRange range;
    if (count==0 || maximum < origin || minimum > origin+delta*count) return range;

    long long lower=static_cast<long long>(floor((minimum-origin)/delta));
    long long upper=static_cast<long long>(ceil((maximum-origin)/delta))-1;
    // A zero-thickness surface exactly on a grid plane receives one cell.
    if (upper<lower) upper=lower;
    lower=max(0LL,min(lower,static_cast<long long>(count)-1));
    upper=max(0LL,min(upper,static_cast<long long>(count)-1));
    if (lower>upper) return range;
    range.minimum=static_cast<size_t>(lower);
    range.maximum=static_cast<size_t>(upper);
    range.valid=true;
    return range;
}

IndexRange xRange(const Triangle& triangle, const VoxelGrid& grid) {
    return coordinateRange(triangle.minimum.x,triangle.maximum.x,
                           grid.originX(),grid.deltaX(),grid.countX());
}
IndexRange yRange(const Triangle& triangle, const VoxelGrid& grid) {
    return coordinateRange(triangle.minimum.y,triangle.maximum.y,
                           grid.originY(),grid.deltaY(),grid.countY());
}
IndexRange zRange(const Triangle& triangle, const VoxelGrid& grid) {
    return coordinateRange(triangle.minimum.z,triangle.maximum.z,
                           grid.originZ(),grid.deltaZ(),grid.countZ());
}

Vec3 cellCentre(const VoxelGrid& grid, size_t i, size_t j, size_t k) {
    return Vec3{grid.originX()+(static_cast<double>(i)+.5)*grid.deltaX(),
                grid.originY()+(static_cast<double>(j)+.5)*grid.deltaY(),
                grid.originZ()+(static_cast<double>(k)+.5)*grid.deltaZ()};
}

void appendSurfaceTriangles(Surface* surface, vector<Triangle>& destination) {
    if (!surface || surface->vertexCount()<3 || surface->getArea()<=0.f) return;

    shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
    RENIndexedFaceSetBuilder builder(indexedFaceSet);
    builder.AddSurface(surface,0);
    const vector<unsigned int>& indices=indexedFaceSet->GetIndices();
    for (size_t index=0;index+2<indices.size();index+=3) {
        const GENPoint a=GENPoint(indexedFaceSet->GetVertex(indices[index]));
        const GENPoint b=GENPoint(indexedFaceSet->GetVertex(indices[index+1]));
        const GENPoint c=GENPoint(indexedFaceSet->GetVertex(indices[index+2]));
        const Triangle triangle=makeTriangle(
                    Vec3{a[0],a[1],a[2]},
                    Vec3{b[0],b[1],b[2]},
                    Vec3{c[0],c[1],c[2]});
        if (!isDegenerate(triangle)) destination.push_back(triangle);
    }
}

struct SceneGeometry {
    struct TreeCanopy {
        vector<Triangle> leaves;
        float attenuation=0.5f;
    };
    vector<vector<Triangle> > buildings;
    vector<Triangle> districtObstacles;
    vector<Triangle> grounds;
    vector<Triangle> floorSupports;
    vector<Vec3> treeGroundSupports;
    vector<TreeCanopy> treeCanopies;
    Vec3 minimum;
    Vec3 maximum;
    bool hasBounds=false;
};

void extendBounds(SceneGeometry& geometry, const Triangle& triangle) {
    if (!geometry.hasBounds) {
        geometry.minimum=triangle.minimum;
        geometry.maximum=triangle.maximum;
        geometry.hasBounds=true;
        return;
    }
    geometry.minimum.x=min(geometry.minimum.x,triangle.minimum.x);
    geometry.minimum.y=min(geometry.minimum.y,triangle.minimum.y);
    geometry.minimum.z=min(geometry.minimum.z,triangle.minimum.z);
    geometry.maximum.x=max(geometry.maximum.x,triangle.maximum.x);
    geometry.maximum.y=max(geometry.maximum.y,triangle.maximum.y);
    geometry.maximum.z=max(geometry.maximum.z,triangle.maximum.z);
}

void extendBounds(SceneGeometry& geometry, const vector<Triangle>& triangles) {
    for (size_t i=0;i<triangles.size();++i) extendBounds(geometry,triangles[i]);
}

SceneGeometry collectGeometry(XmlScene& scene, const VoxelizationOptions& options) {
    SceneGeometry geometry;
    District* district=scene.getDistrict();
    if (!district) throw runtime_error("Voxelization: the scene has no District.");

    for (size_t buildingIndex=0;buildingIndex<district->getnBuildings();++buildingIndex) {
        Building* building=district->getBuilding(buildingIndex);
        vector<Triangle> triangles;
        for (size_t zoneIndex=0;zoneIndex<building->getnZones();++zoneIndex) {
            Zone* zone=building->getZone(zoneIndex);
            const vector<Surface*> surfaces=zone->getAllSurfaces();
            for (size_t surfaceIndex=0;surfaceIndex<surfaces.size();++surfaceIndex)
                appendSurfaceTriangles(surfaces[surfaceIndex],triangles);
            for (size_t floorIndex=0;floorIndex<zone->getnFloors();++floorIndex)
                appendSurfaceTriangles(zone->getFloor(floorIndex),geometry.floorSupports);
        }
        if (!triangles.empty()) {
            extendBounds(geometry,triangles);
            geometry.buildings.push_back(move(triangles));
        }
    }

    if (options.includeDistrictSurfacesAsSolids) {
        for (size_t surfaceIndex=0;surfaceIndex<district->getnSurfaces();++surfaceIndex)
            appendSurfaceTriangles(district->getSurface(surfaceIndex),geometry.districtObstacles);
        extendBounds(geometry,geometry.districtObstacles);

        const vector<Surface*>& nonSimulatedFloors=
                district->getNonSimulatedBuildingFloors();
        for (size_t floorIndex=0;floorIndex<nonSimulatedFloors.size();++floorIndex)
            appendSurfaceTriangles(nonSimulatedFloors[floorIndex],
                                   geometry.floorSupports);
    }

    for (forward_list<Ground*>::iterator ground=district->getGrounds()->begin();
         ground!=district->getGrounds()->end();++ground)
        appendSurfaceTriangles(*ground,geometry.grounds);
    extendBounds(geometry,geometry.grounds);

    for (size_t treeIndex=0;treeIndex<district->getnTrees();++treeIndex) {
        Tree* tree=district->getTree(treeIndex);
        vector<Triangle> leaves;
        // Include explicit leaves and generated layers; Trunc is never material.
        for (size_t leafIndex=0;leafIndex<tree->getLeaves()->size();++leafIndex)
            appendSurfaceTriangles(tree->getLeaves()->at(leafIndex),leaves);
        for (size_t leafIndex=0;leafIndex<tree->getSubLeaves()->size();++leafIndex)
            appendSurfaceTriangles(tree->getSubLeaves()->at(leafIndex),leaves);
        if (!leaves.empty()) {
            extendBounds(geometry,leaves);
            SceneGeometry::TreeCanopy canopy;
            canopy.leaves=move(leaves);
            canopy.attenuation=tree->getWindAttenuation();
            geometry.treeCanopies.push_back(move(canopy));
        }

        // Trunks remain absent from the material grid, but their lowest
        // vertices are reliable terrain-height evidence in tree-only regions.
        double baseZ=numeric_limits<double>::infinity();
        double baseX=0.,baseY=0.;
        size_t baseVertexCount=0;
        for (size_t truncIndex=0;truncIndex<tree->getTrunc()->size();++truncIndex) {
            const vector<GENPoint>& vertices=
                    *tree->getTrunc()->at(truncIndex)->getVertices();
            for (size_t vertexIndex=0;vertexIndex<vertices.size();++vertexIndex) {
                const double z=vertices[vertexIndex][2];
                if (z<baseZ-1e-6) {
                    baseZ=z;
                    baseX=vertices[vertexIndex][0];
                    baseY=vertices[vertexIndex][1];
                    baseVertexCount=1;
                }
                else if (fabs(z-baseZ)<=1e-6) {
                    baseX+=vertices[vertexIndex][0];
                    baseY+=vertices[vertexIndex][1];
                    ++baseVertexCount;
                }
            }
        }
        if (baseVertexCount>0) {
            geometry.treeGroundSupports.push_back(
                        Vec3{baseX/baseVertexCount,baseY/baseVertexCount,baseZ});
        }
    }

    if (!geometry.hasBounds)
        throw runtime_error("Voxelization: no usable building, ground, obstacle or leaf triangles.");
    for (size_t supportIndex=0;
         supportIndex<geometry.treeGroundSupports.size();++supportIndex)
        geometry.minimum.z=min(geometry.minimum.z,
                               geometry.treeGroundSupports[supportIndex].z);
    return geometry;
}

int effectiveThreadCount(const VoxelizationOptions& options) {
#ifdef _OPENMP
    return options.threads>0 ? options.threads : omp_get_max_threads();
#else
    (void)options;
    return 1;
#endif
}

size_t checkedCellCount(size_t nx, size_t ny, size_t nz,
                             size_t maximum) {
    if (nx==0 || ny==0 || nz==0) throw runtime_error("Voxelization: empty grid.");
    if (ny>numeric_limits<size_t>::max()/nx ||
        nz>numeric_limits<size_t>::max()/(nx*ny))
        throw runtime_error("Voxelization: grid dimensions overflow size_t.");
    const size_t count=nx*ny*nz;
    if (maximum>0 && count>maximum)
        throw runtime_error("Voxelization: requested grid exceeds maximumCellCount.");
    if (count>numeric_limits<uint32_t>::max())
        throw runtime_error("Voxelization: V1 exterior fill supports at most 2^32-1 cells.");
    return count;
}

struct LocalMask {
    size_t i0=0,j0=0,k0=0;
    size_t nx=0,ny=0,nz=0;
    vector<uint8_t> values;
    size_t index(size_t i,size_t j,size_t k) const {
        return (k*ny+j)*nx+i;
    }
};

LocalMask rasterizeBuilding(const vector<Triangle>& triangles,
                            const VoxelGrid& grid,
                            unsigned int closingPasses) {
    LocalMask mask;
    if (triangles.empty()) return mask;
    size_t minI=grid.countX()-1,minJ=grid.countY()-1,minK=grid.countZ()-1;
    size_t maxI=0,maxJ=0,maxK=0;
    bool valid=false;
    for (size_t triangleIndex=0;triangleIndex<triangles.size();++triangleIndex) {
        const IndexRange xr=xRange(triangles[triangleIndex],grid);
        const IndexRange yr=yRange(triangles[triangleIndex],grid);
        const IndexRange zr=zRange(triangles[triangleIndex],grid);
        if (!xr.valid || !yr.valid || !zr.valid) continue;
        minI=min(minI,xr.minimum); maxI=max(maxI,xr.maximum);
        minJ=min(minJ,yr.minimum); maxJ=max(maxJ,yr.maximum);
        minK=min(minK,zr.minimum); maxK=max(maxK,zr.maximum);
        valid=true;
    }
    if (!valid) return mask;

    const size_t padding=static_cast<size_t>(closingPasses)+1;
    mask.i0=minI>padding?minI-padding:0;
    mask.j0=minJ>padding?minJ-padding:0;
    mask.k0=minK>padding?minK-padding:0;
    maxI=min(grid.countX()-1,maxI+padding);
    maxJ=min(grid.countY()-1,maxJ+padding);
    maxK=min(grid.countZ()-1,maxK+padding);
    mask.nx=maxI-mask.i0+1;
    mask.ny=maxJ-mask.j0+1;
    mask.nz=maxK-mask.k0+1;
    mask.values.assign(mask.nx*mask.ny*mask.nz,0);

    const Vec3 half{.5*grid.deltaX(),.5*grid.deltaY(),.5*grid.deltaZ()};
    for (size_t triangleIndex=0;triangleIndex<triangles.size();++triangleIndex) {
        const Triangle& triangle=triangles[triangleIndex];
        const IndexRange xr=xRange(triangle,grid);
        const IndexRange yr=yRange(triangle,grid);
        const IndexRange zr=zRange(triangle,grid);
        if (!xr.valid || !yr.valid || !zr.valid) continue;
        for (size_t k=zr.minimum;k<=zr.maximum;++k)
            for (size_t j=yr.minimum;j<=yr.maximum;++j)
                for (size_t i=xr.minimum;i<=xr.maximum;++i)
                    if (triangleIntersectsBox(triangle,cellCentre(grid,i,j,k),half))
                        mask.values[mask.index(i-mask.i0,j-mask.j0,k-mask.k0)]=1;
    }

    // Six-neighbour closing is consistent with the later finite-volume stencil.
    for (unsigned int pass=0;pass<closingPasses;++pass) {
        vector<uint8_t> dilated(mask.values.size(),0);
        for (size_t k=0;k<mask.nz;++k)
            for (size_t j=0;j<mask.ny;++j)
                for (size_t i=0;i<mask.nx;++i) {
                    const size_t index=mask.index(i,j,k);
                    if (mask.values[index] ||
                        (i>0 && mask.values[mask.index(i-1,j,k)]) ||
                        (i+1<mask.nx && mask.values[mask.index(i+1,j,k)]) ||
                        (j>0 && mask.values[mask.index(i,j-1,k)]) ||
                        (j+1<mask.ny && mask.values[mask.index(i,j+1,k)]) ||
                        (k>0 && mask.values[mask.index(i,j,k-1)]) ||
                        (k+1<mask.nz && mask.values[mask.index(i,j,k+1)]))
                        dilated[index]=1;
                }

        // A closing operation must be extensive: never erase an originally
        // intersected shell cell, including shells touching the scene AABB.
        vector<uint8_t> closed(mask.values);
        for (size_t k=1;k+1<mask.nz;++k)
            for (size_t j=1;j+1<mask.ny;++j)
                for (size_t i=1;i+1<mask.nx;++i) {
                    const size_t index=mask.index(i,j,k);
                    if (dilated[index] &&
                        dilated[mask.index(i-1,j,k)] && dilated[mask.index(i+1,j,k)] &&
                        dilated[mask.index(i,j-1,k)] && dilated[mask.index(i,j+1,k)] &&
                        dilated[mask.index(i,j,k-1)] && dilated[mask.index(i,j,k+1)])
                        closed[index]=1;
                }
        mask.values.swap(closed);
    }
    return mask;
}

void mergeBuildings(const SceneGeometry& geometry, VoxelGrid& grid,
                    const VoxelizationOptions& options, int threadCount) {
    vector<uint8_t>& materials=grid.materials();
    const long long buildingCount=static_cast<long long>(geometry.buildings.size());

    #pragma omp parallel for schedule(dynamic) num_threads(threadCount) if(buildingCount>1)
    for (long long buildingIndex=0;buildingIndex<buildingCount;++buildingIndex) {
        const LocalMask mask=rasterizeBuilding(
                    geometry.buildings[static_cast<size_t>(buildingIndex)],
                    grid,options.buildingGapClosingPasses);
        for (size_t k=0;k<mask.nz;++k)
            for (size_t j=0;j<mask.ny;++j)
                for (size_t i=0;i<mask.nx;++i)
                    if (mask.values[mask.index(i,j,k)]) {
                        const size_t global=grid.index(mask.i0+i,mask.j0+j,mask.k0+k);
                        #pragma omp atomic write
                        materials[global]=BUILDING_VALUE;
                    }
    }

    // Generic District surfaces have no reliable per-building grouping.
    const Vec3 half{.5*grid.deltaX(),.5*grid.deltaY(),.5*grid.deltaZ()};
    const long long obstacleCount=static_cast<long long>(geometry.districtObstacles.size());
    #pragma omp parallel for schedule(dynamic) num_threads(threadCount) if(obstacleCount>32)
    for (long long obstacleIndex=0;obstacleIndex<obstacleCount;++obstacleIndex) {
        const Triangle& triangle=geometry.districtObstacles[static_cast<size_t>(obstacleIndex)];
        const IndexRange xr=xRange(triangle,grid);
        const IndexRange yr=yRange(triangle,grid);
        const IndexRange zr=zRange(triangle,grid);
        if (!xr.valid || !yr.valid || !zr.valid) continue;
        for (size_t k=zr.minimum;k<=zr.maximum;++k)
            for (size_t j=yr.minimum;j<=yr.maximum;++j)
                for (size_t i=xr.minimum;i<=xr.maximum;++i)
                    if (triangleIntersectsBox(triangle,cellCentre(grid,i,j,k),half)) {
                        const size_t index=grid.index(i,j,k);
                        #pragma omp atomic write
                        materials[index]=BUILDING_VALUE;
                    }
    }
}

bool projectedPointHeight(const Triangle& triangle, double x, double y, double& z) {
    const double x0=triangle.p[0].x,y0=triangle.p[0].y;
    const double x1=triangle.p[1].x,y1=triangle.p[1].y;
    const double x2=triangle.p[2].x,y2=triangle.p[2].y;
    const double denominator=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    if (fabs(denominator)<1e-18) return false;
    const double a=((y1-y2)*(x-x2)+(x2-x1)*(y-y2))/denominator;
    const double b=((y2-y0)*(x-x2)+(x0-x2)*(y-y2))/denominator;
    const double c=1.-a-b;
    const double tolerance=1e-10;
    if (a < -tolerance || b < -tolerance || c < -tolerance) return false;
    z=a*triangle.p[0].z+b*triangle.p[1].z+c*triangle.p[2].z;
    return true;
}

bool fillGroundHeightHoles(vector<float>& heights,
                           size_t nx,size_t ny,int threadCount) {
    vector<uint8_t> rowHasValue(ny,0);
    const long long rowCount=static_cast<long long>(ny);

    #pragma omp parallel for schedule(static) num_threads(threadCount) if(rowCount>16)
    for (long long rowValue=0;rowValue<rowCount;++rowValue) {
        const size_t row=static_cast<size_t>(rowValue);
        size_t first=0;
        while (first<nx && !isfinite(heights[row*nx+first])) ++first;
        if (first==nx) continue;
        rowHasValue[row]=1;

        for (size_t i=0;i<first;++i) heights[row*nx+i]=heights[row*nx+first];
        size_t previous=first;
        for (size_t next=first+1;next<nx;++next) {
            if (!isfinite(heights[row*nx+next])) continue;
            const float left=heights[row*nx+previous];
            const float right=heights[row*nx+next];
            for (size_t i=previous+1;i<next;++i) {
                const double fraction=static_cast<double>(i-previous)/static_cast<double>(next-previous);
                heights[row*nx+i]=static_cast<float>(left+(right-left)*fraction);
            }
            previous=next;
        }
        for (size_t i=previous+1;i<nx;++i) heights[row*nx+i]=heights[row*nx+previous];
    }

    size_t firstRow=0;
    while (firstRow<ny && !rowHasValue[firstRow]) ++firstRow;
    if (firstRow==ny) return false;

    for (size_t row=0;row<firstRow;++row)
        copy(heights.begin()+firstRow*nx,heights.begin()+(firstRow+1)*nx,
                  heights.begin()+row*nx);

    size_t previousRow=firstRow;
    for (size_t nextRow=firstRow+1;nextRow<ny;++nextRow) {
        if (!rowHasValue[nextRow]) continue;
        for (size_t row=previousRow+1;row<nextRow;++row) {
            const double fraction=static_cast<double>(row-previousRow)/
                                  static_cast<double>(nextRow-previousRow);
            for (size_t i=0;i<nx;++i) {
                const float lower=heights[previousRow*nx+i];
                const float upper=heights[nextRow*nx+i];
                heights[row*nx+i]=static_cast<float>(lower+(upper-lower)*fraction);
            }
        }
        previousRow=nextRow;
    }
    for (size_t row=previousRow+1;row<ny;++row)
        copy(heights.begin()+previousRow*nx,heights.begin()+(previousRow+1)*nx,
                  heights.begin()+row*nx);
    return true;
}

vector<vector<uint32_t> > makeHorizontalBins(
        const vector<Triangle>& triangles,
        const VoxelGrid& grid,
        size_t tileSize,
        size_t tileCountX,
        size_t tileCountY) {
    vector<vector<uint32_t> > bins(tileCountX*tileCountY);
    for (size_t triangleIndex=0;triangleIndex<triangles.size();++triangleIndex) {
        const IndexRange xr=xRange(triangles[triangleIndex],grid);
        const IndexRange yr=yRange(triangles[triangleIndex],grid);
        if (!xr.valid || !yr.valid) continue;
        const size_t tileMinX=xr.minimum/tileSize,tileMaxX=xr.maximum/tileSize;
        const size_t tileMinY=yr.minimum/tileSize,tileMaxY=yr.maximum/tileSize;
        for (size_t tileY=tileMinY;tileY<=tileMaxY;++tileY)
            for (size_t tileX=tileMinX;tileX<=tileMaxX;++tileX)
                bins[tileY*tileCountX+tileX].push_back(
                            static_cast<uint32_t>(triangleIndex));
    }
    return bins;
}

struct GroundReconstructionStats {
    size_t explicitGroundColumns=0;
    size_t buildingFloorColumns=0;
    size_t treeBaseColumns=0;
    bool usedGlobalMinimum=false;
};

GroundReconstructionStats voxelizeGround(const SceneGeometry& geometry,
                                         VoxelGrid& grid,
                                         const VoxelizationOptions& options,
                                         int threadCount) {
    GroundReconstructionStats stats;
    const size_t nx=grid.countX(),ny=grid.countY();
    const size_t tileSize=max<size_t>(1,options.groundTileSize);
    const size_t tileCountX=(nx+tileSize-1)/tileSize;
    const size_t tileCountY=(ny+tileSize-1)/tileSize;
    const vector<vector<uint32_t> > groundBins=makeHorizontalBins(
                geometry.grounds,grid,tileSize,tileCountX,tileCountY);
    const vector<vector<uint32_t> > floorBins=makeHorizontalBins(
                geometry.floorSupports,grid,tileSize,tileCountX,tileCountY);

    vector<float>& heights=grid.groundHeights();
    vector<uint8_t> source(nx*ny,0);
    const Vec3 columnHalf{.5*grid.deltaX(),.5*grid.deltaY(),
                          .5*(grid.maximumZ()-grid.originZ())+grid.deltaZ()};
    const double columnCentreZ=.5*(grid.maximumZ()+grid.originZ());
    const long long columnCount=static_cast<long long>(nx*ny);
    long long explicitGroundColumns=0;
    long long buildingFloorColumns=0;

    #pragma omp parallel for schedule(static) num_threads(threadCount) if(columnCount>1024) \
        reduction(+:explicitGroundColumns,buildingFloorColumns)
    for (long long columnValue=0;columnValue<columnCount;++columnValue) {
        const size_t column=static_cast<size_t>(columnValue);
        const size_t i=column%nx,j=column/nx;
        const double x=grid.originX()+(static_cast<double>(i)+.5)*grid.deltaX();
        const double y=grid.originY()+(static_cast<double>(j)+.5)*grid.deltaY();
        const size_t bin=(j/tileSize)*tileCountX+(i/tileSize);
        const vector<uint32_t>& candidates=groundBins[bin];
        double maximumHeight=-numeric_limits<double>::infinity();
        for (size_t candidate=0;candidate<candidates.size();++candidate) {
            const Triangle& triangle=geometry.grounds[candidates[candidate]];
            double z=0.;
            if (projectedPointHeight(triangle,x,y,z)) {
                maximumHeight=max(maximumHeight,z);
            }
            else if (triangleIntersectsBox(
                         triangle,Vec3{x,y,columnCentreZ},columnHalf)) {
                // The surface crosses this column but misses its centre.  Using
                // its highest vertex is the conservative LEGO approximation.
                maximumHeight=max(maximumHeight,triangle.maximum.z);
            }
        }
        if (isfinite(maximumHeight)) {
            heights[column]=static_cast<float>(maximumHeight);
            source[column]=1;
            ++explicitGroundColumns;
            continue;
        }

        // A lowest building floor is a local terrain observation when the
        // source scene provides no explicit ground below that footprint.
        const vector<uint32_t>& floorCandidates=floorBins[bin];
        double minimumHeight=numeric_limits<double>::infinity();
        for (size_t candidate=0;candidate<floorCandidates.size();++candidate) {
            const Triangle& triangle=
                    geometry.floorSupports[floorCandidates[candidate]];
            double z=0.;
            if (projectedPointHeight(triangle,x,y,z)) {
                minimumHeight=min(minimumHeight,z);
            }
            else if (triangleIntersectsBox(
                         triangle,Vec3{x,y,columnCentreZ},columnHalf)) {
                minimumHeight=min(minimumHeight,triangle.minimum.z);
            }
        }
        if (isfinite(minimumHeight)) {
            heights[column]=static_cast<float>(minimumHeight);
            source[column]=2;
            ++buildingFloorColumns;
        }
    }

    stats.explicitGroundColumns=static_cast<size_t>(explicitGroundColumns);
    stats.buildingFloorColumns=static_cast<size_t>(buildingFloorColumns);

    // Tree trunks are not solids, but their base is useful as a sparse ground
    // elevation observation when neither ground nor a floor covers the column.
    for (size_t supportIndex=0;
         supportIndex<geometry.treeGroundSupports.size();++supportIndex) {
        const Vec3& support=geometry.treeGroundSupports[supportIndex];
        long long i=static_cast<long long>(
                    floor((support.x-grid.originX())/grid.deltaX()));
        long long j=static_cast<long long>(
                    floor((support.y-grid.originY())/grid.deltaY()));
        i=max(0LL,min(i,static_cast<long long>(nx)-1));
        j=max(0LL,min(j,static_cast<long long>(ny)-1));
        const size_t column=static_cast<size_t>(j)*nx+static_cast<size_t>(i);
        if (source[column]==0) {
            heights[column]=static_cast<float>(support.z);
            source[column]=3;
            ++stats.treeBaseColumns;
        }
        else if (source[column]==3) {
            heights[column]=max(heights[column],static_cast<float>(support.z));
        }
    }

    if (!fillGroundHeightHoles(heights,nx,ny,threadCount)) {
        fill(heights.begin(),heights.end(),static_cast<float>(geometry.minimum.z));
        stats.usedGlobalMinimum=true;
    }

    vector<uint8_t>& materials=grid.materials();
    #pragma omp parallel for schedule(static) num_threads(threadCount) if(columnCount>1024)
    for (long long columnValue=0;columnValue<columnCount;++columnValue) {
        const size_t column=static_cast<size_t>(columnValue);
        const size_t i=column%nx,j=column/nx;
        long long surfaceK=static_cast<long long>(
                    floor((static_cast<double>(heights[column])-grid.originZ())/
                               grid.deltaZ()));
        surfaceK=max(0LL,min(surfaceK,static_cast<long long>(grid.countZ())-1));
        for (long long k=0;k<=surfaceK;++k)
            materials[grid.index(i,j,static_cast<size_t>(k))]=GROUND_VALUE;
    }
    return stats;
}

void fillEnclosedBuildingCells(VoxelGrid& grid, bool hasBuildingShell, int threadCount) {
    if (!hasBuildingShell) return;
    vector<uint8_t>& materials=grid.materials();
    vector<uint32_t> stack;
    const size_t boundaryEstimate=
            2*(grid.countX()*grid.countY()+grid.countX()*grid.countZ()+grid.countY()*grid.countZ());
    stack.reserve(min<size_t>(boundaryEstimate,grid.cellCount()));

    const auto seed=[&](size_t i,size_t j,size_t k) {
        const size_t index=grid.index(i,j,k);
        if (materials[index]==AIR_VALUE) {
            materials[index]=OUTSIDE_AIR_VALUE;
            stack.push_back(static_cast<uint32_t>(index));
        }
    };
    for (size_t k=0;k<grid.countZ();++k)
        for (size_t j=0;j<grid.countY();++j) {
            seed(0,j,k); seed(grid.countX()-1,j,k);
        }
    for (size_t k=0;k<grid.countZ();++k)
        for (size_t i=0;i<grid.countX();++i) {
            seed(i,0,k); seed(i,grid.countY()-1,k);
        }
    for (size_t j=0;j<grid.countY();++j)
        for (size_t i=0;i<grid.countX();++i) {
            seed(i,j,0); seed(i,j,grid.countZ()-1);
        }

    const auto visit=[&](size_t index) {
        if (materials[index]==AIR_VALUE) {
            materials[index]=OUTSIDE_AIR_VALUE;
            stack.push_back(static_cast<uint32_t>(index));
        }
    };
    const size_t slice=grid.countX()*grid.countY();
    while (!stack.empty()) {
        const size_t index=stack.back();
        stack.pop_back();
        const size_t i=index%grid.countX();
        const size_t plane=index/grid.countX();
        const size_t j=plane%grid.countY();
        const size_t k=plane/grid.countY();
        if (i>0) visit(index-1);
        if (i+1<grid.countX()) visit(index+1);
        if (j>0) visit(index-grid.countX());
        if (j+1<grid.countY()) visit(index+grid.countX());
        if (k>0) visit(index-slice);
        if (k+1<grid.countZ()) visit(index+slice);
    }

    const long long cellCount=static_cast<long long>(materials.size());
    #pragma omp parallel for schedule(static) num_threads(threadCount) if(cellCount>100000)
    for (long long index=0;index<cellCount;++index) {
        if (materials[static_cast<size_t>(index)]==AIR_VALUE)
            materials[static_cast<size_t>(index)]=BUILDING_VALUE;
        else if (materials[static_cast<size_t>(index)]==OUTSIDE_AIR_VALUE)
            materials[static_cast<size_t>(index)]=AIR_VALUE;
    }
}

void voxelizeTrees(const vector<SceneGeometry::TreeCanopy>& trees,
                   VoxelGrid& grid,int threadCount) {
    const size_t columnCount=grid.countX()*grid.countY();
    vector<float> minimum(columnCount,numeric_limits<float>::infinity());
    vector<float> maximum(columnCount,-numeric_limits<float>::infinity());
    vector<uint32_t> touched;
    const Vec3 columnHalf{.5*grid.deltaX(),.5*grid.deltaY(),
                          .5*(grid.maximumZ()-grid.originZ())+grid.deltaZ()};
    const double columnCentreZ=.5*(grid.maximumZ()+grid.originZ());
    vector<uint8_t>& materials=grid.materials();
    vector<pair<uint32_t,float> > attenuationAssignments;

    for (size_t treeIndex=0;treeIndex<trees.size();++treeIndex) {
        touched.clear();
        const vector<Triangle>& leaves=trees[treeIndex].leaves;
        for (size_t triangleIndex=0;triangleIndex<leaves.size();++triangleIndex) {
            const Triangle& triangle=leaves[triangleIndex];
            const IndexRange xr=xRange(triangle,grid);
            const IndexRange yr=yRange(triangle,grid);
            if (!xr.valid || !yr.valid) continue;
            for (size_t j=yr.minimum;j<=yr.maximum;++j)
                for (size_t i=xr.minimum;i<=xr.maximum;++i) {
                    const Vec3 centre{
                        grid.originX()+(static_cast<double>(i)+.5)*grid.deltaX(),
                        grid.originY()+(static_cast<double>(j)+.5)*grid.deltaY(),
                        columnCentreZ};
                    if (!triangleIntersectsBox(triangle,centre,columnHalf)) continue;
                    const size_t column=j*grid.countX()+i;
                    if (!isfinite(minimum[column]))
                        touched.push_back(static_cast<uint32_t>(column));
                    minimum[column]=min(minimum[column],static_cast<float>(triangle.minimum.z));
                    maximum[column]=max(maximum[column],static_cast<float>(triangle.maximum.z));
                }
        }

        const long long touchedCount=static_cast<long long>(touched.size());
        #pragma omp parallel for schedule(static) num_threads(threadCount) if(touchedCount>256)
        for (long long touchedIndex=0;touchedIndex<touchedCount;++touchedIndex) {
            const size_t column=touched[static_cast<size_t>(touchedIndex)];
            const size_t i=column%grid.countX(),j=column/grid.countX();
            long long kMinimum=static_cast<long long>(
                        floor((minimum[column]-grid.originZ())/grid.deltaZ()));
            long long kMaximum=static_cast<long long>(
                        floor((maximum[column]-grid.originZ())/grid.deltaZ()));
            kMinimum=max(0LL,min(kMinimum,static_cast<long long>(grid.countZ())-1));
            kMaximum=max(0LL,min(kMaximum,static_cast<long long>(grid.countZ())-1));
            for (long long k=kMinimum;k<=kMaximum;++k) {
                const size_t index=grid.index(i,j,static_cast<size_t>(k));
                if (materials[index]==AIR_VALUE) materials[index]=TREE_VALUE;
            }
        }

        // Preserve only sparse TREE-cell assignments. Overlapping canopies
        // are merged later with the strongest (smallest) attenuation.
        for (size_t touchedIndex=0;touchedIndex<touched.size();++touchedIndex) {
            const size_t column=touched[touchedIndex];
            const size_t i=column%grid.countX(),j=column/grid.countX();
            long long kMinimum=static_cast<long long>(
                        floor((minimum[column]-grid.originZ())/grid.deltaZ()));
            long long kMaximum=static_cast<long long>(
                        floor((maximum[column]-grid.originZ())/grid.deltaZ()));
            kMinimum=max(0LL,min(kMinimum,static_cast<long long>(grid.countZ())-1));
            kMaximum=max(0LL,min(kMaximum,static_cast<long long>(grid.countZ())-1));
            for (long long k=kMinimum;k<=kMaximum;++k) {
                const size_t cell=grid.index(i,j,static_cast<size_t>(k));
                if (materials[cell]==TREE_VALUE)
                    attenuationAssignments.push_back(make_pair(
                        static_cast<uint32_t>(cell),trees[treeIndex].attenuation));
            }
        }

        for (size_t touchedIndex=0;touchedIndex<touched.size();++touchedIndex) {
            const size_t column=touched[touchedIndex];
            minimum[column]=numeric_limits<float>::infinity();
            maximum[column]=-numeric_limits<float>::infinity();
        }
    }
    grid.setTreeWindAttenuations(move(attenuationAssignments));
}

} // namespace

VoxelGrid::VoxelGrid() {}

VoxelGrid::VoxelGrid(double originX, double originY, double originZ,
                     double deltaX, double deltaY, double deltaZ,
                     size_t countX, size_t countY, size_t countZ)
    :originX_(originX),originY_(originY),originZ_(originZ),
      deltaX_(deltaX),deltaY_(deltaY),deltaZ_(deltaZ),
      countX_(countX),countY_(countY),countZ_(countZ),
      materials_(countX*countY*countZ,AIR_VALUE),
      groundHeights_(countX*countY,numeric_limits<float>::quiet_NaN()) {}

void VoxelGrid::setLogLawWindProfile(double referenceSpeed,
                                     double meteorologicalDirectionDegrees,
                                     double roughnessLength,
                                     double referenceHeight) {
    if (!isfinite(referenceSpeed) || referenceSpeed<0.)
        throw runtime_error("Voxelization: climate reference wind speed must be finite and non-negative.");
    if (!isfinite(meteorologicalDirectionDegrees))
        throw runtime_error("Voxelization: climate wind direction must be finite.");
    if (!isfinite(roughnessLength) || !(roughnessLength>0.)
            || !isfinite(referenceHeight) || !(referenceHeight>roughnessLength))
        throw runtime_error("Voxelization: roughness must be finite, positive and below the wind reference height.");

    double direction=fmod(meteorologicalDirectionDegrees,360.);
    if (direction<0.) direction+=360.;
    const double radians=direction*acos(-1.)/180.;
    const double denominator=log1p(referenceHeight/roughnessLength);

    windSpeedByLevel_.resize(countZ_);
    windVelocityByLevel_.resize(countZ_);
    for (size_t k=0;k<countZ_;++k) {
        // The unattenuated V1 profile is horizontally uniform: height is
        // measured from the lower domain plane at each voxel centre.
        const double height=(static_cast<double>(k)+.5)*deltaZ_;
        const double speed=referenceSpeed
                *log1p(height/roughnessLength)/denominator;
        windSpeedByLevel_[k]=static_cast<float>(speed);
        // Climate DD is meteorological: clockwise from North and indicating
        // where the wind comes from. CitySim x is East and y is North.
        windVelocityByLevel_[k]={{
                static_cast<float>(-speed*sin(radians)),
                static_cast<float>(-speed*cos(radians)),
                0.f
        }};
    }

    referenceWindSpeed_=referenceSpeed;
    referenceWindDirection_=direction;
    roughnessLength_=roughnessLength;
    referenceWindHeight_=referenceHeight;
}

void VoxelGrid::setTreeWindAttenuations(
        vector<pair<uint32_t,float> > assignments) {
    sort(assignments.begin(),assignments.end(),
         [](const pair<uint32_t,float>& first,
            const pair<uint32_t,float>& second) {
        if (first.first!=second.first) return first.first<second.first;
        return first.second<second.second;
    });

    treeWindCellIndices_.clear();
    treeWindAttenuations_.clear();
    treeWindCellIndices_.reserve(assignments.size());
    treeWindAttenuations_.reserve(assignments.size());
    for (size_t assignment=0;assignment<assignments.size();++assignment) {
        const uint32_t cell=assignments[assignment].first;
        const float attenuation=assignments[assignment].second;
        if (cell>=materials_.size() || materials_[cell]!=TREE_VALUE)
            throw runtime_error("Voxelization: tree wind assignment targets a non-TREE cell.");
        if (!isfinite(attenuation) || !(attenuation>0.f)
                || !(attenuation<1.f))
            throw runtime_error("Voxelization: Tree attenuation must be finite and strictly between zero and one.");

        if (!treeWindCellIndices_.empty()
                && treeWindCellIndices_.back()==cell) {
            treeWindAttenuations_.back()=
                    min(treeWindAttenuations_.back(),attenuation);
        }
        else {
            treeWindCellIndices_.push_back(cell);
            treeWindAttenuations_.push_back(attenuation);
        }
    }

    const size_t treeCells=static_cast<size_t>(
                count(materials_.begin(),materials_.end(),TREE_VALUE));
    if (treeWindCellIndices_.size()!=treeCells)
        throw runtime_error("Voxelization: not every TREE cell received a wind attenuation.");
}

VoxelCounts VoxelGrid::counts() const {
    VoxelCounts result;
    for (size_t index=0;index<materials_.size();++index) {
        switch (static_cast<VoxelMaterial>(materials_[index])) {
        case VoxelMaterial::AIR: ++result.air; break;
        case VoxelMaterial::BUILDING: ++result.building; break;
        case VoxelMaterial::GROUND: ++result.ground; break;
        case VoxelMaterial::TREE: ++result.tree; break;
        }
    }
    return result;
}

void VoxelGrid::writeLegacyVtk(const string& fileName) const {
    ofstream output(fileName.c_str(),ios::out|ios::binary|ios::trunc);
    if (!output.good()) throw runtime_error("Voxelization: cannot create VTK file "+fileName);
    output << "# vtk DataFile Version 3.0\n"
           << "CitySim voxel materials: 0=AIR 1=BUILDING 2=GROUND 3=TREE\n"
           << "BINARY\n"
           << "DATASET STRUCTURED_POINTS\n"
           << "DIMENSIONS " << countX_+1 << " " << countY_+1 << " " << countZ_+1 << "\n"
           << setprecision(17)
           << "ORIGIN " << originX_ << " " << originY_ << " " << originZ_ << "\n"
           << "SPACING " << deltaX_ << " " << deltaY_ << " " << deltaZ_ << "\n"
           << "CELL_DATA " << materials_.size() << "\n"
           << "SCALARS material unsigned_char 1\n"
           << "LOOKUP_TABLE default\n";
    if (!materials_.empty())
        output.write(reinterpret_cast<const char*>(&materials_[0]),
                     static_cast<streamsize>(materials_.size()));
    output << "\n";

    if (hasWindProfile()) {
        // Legacy binary VTK requires multi-byte numbers in big-endian order.
        const auto storeFloat = [](vector<char>& buffer,size_t offset,float value) {
            uint32_t bits=0;
            memcpy(&bits,&value,sizeof(bits));
            buffer[offset]=static_cast<char>((bits>>24)&0xffu);
            buffer[offset+1]=static_cast<char>((bits>>16)&0xffu);
            buffer[offset+2]=static_cast<char>((bits>>8)&0xffu);
            buffer[offset+3]=static_cast<char>(bits&0xffu);
        };
        const auto attenuationAt =
                [this](size_t cell,size_t& treeCursor)->float {
            if (materials_[cell]!=TREE_VALUE) return 1.f;
            while (treeCursor<treeWindCellIndices_.size()
                   && treeWindCellIndices_[treeCursor]<cell)
                ++treeCursor;
            if (treeCursor>=treeWindCellIndices_.size()
                    || treeWindCellIndices_[treeCursor]!=cell)
                throw runtime_error(
                    "Voxelization: TREE cell has no wind attenuation.");
            return treeWindAttenuations_[treeCursor++];
        };

        output << "SCALARS wind_speed float 1\n"
               << "LOOKUP_TABLE default\n";
        vector<char> scalarRow(countX_*sizeof(float));
        size_t treeCursor=0;
        for (size_t k=0;k<countZ_;++k) {
            for (size_t j=0;j<countY_;++j) {
                for (size_t i=0;i<countX_;++i) {
                    const size_t cell=index(i,j,k);
                    const float attenuation=attenuationAt(cell,treeCursor);
                    storeFloat(scalarRow,i*sizeof(float),
                               attenuation*windSpeedByLevel_[k]);
                }
                output.write(&scalarRow[0],static_cast<streamsize>(scalarRow.size()));
            }
        }
        output << "\nVECTORS wind_velocity float\n";
        vector<char> vectorRow(countX_*3*sizeof(float));
        treeCursor=0;
        for (size_t k=0;k<countZ_;++k) {
            for (size_t j=0;j<countY_;++j) {
                for (size_t i=0;i<countX_;++i) {
                    const size_t cell=index(i,j,k);
                    const float attenuation=attenuationAt(cell,treeCursor);
                    for (size_t component=0;component<3;++component)
                        storeFloat(vectorRow,(i*3+component)*sizeof(float),
                                   attenuation
                                   *windVelocityByLevel_[k][component]);
                }
                output.write(&vectorRow[0],static_cast<streamsize>(vectorRow.size()));
            }
        }
        output << "\n";
    }
    if (!output.good()) throw runtime_error("Voxelization: error writing VTK file "+fileName);
}

VoxelGrid SceneVoxelizer::voxelize(XmlScene& scene,
                                   const VoxelizationOptions& options,
                                   ostream* requestedLog) {
    const chrono::steady_clock::time_point start=chrono::steady_clock::now();
    ostream& log=requestedLog?*requestedLog:cout;
    if (!(options.deltaX>0.) || !(options.deltaY>0.) || !(options.deltaZ>0.) ||
        !isfinite(options.deltaX) || !isfinite(options.deltaY) ||
        !isfinite(options.deltaZ))
        throw runtime_error("Voxelization: deltaX, deltaY and deltaZ must be finite and positive.");
    if (options.groundTileSize==0)
        throw runtime_error("Voxelization: groundTileSize must be greater than zero.");
    if (options.threads<0)
        throw runtime_error("Voxelization: threads cannot be negative.");
    if (!isfinite(options.roughnessLength) || !(options.roughnessLength>0.)
            || !(options.roughnessLength<10.))
        throw runtime_error("Voxelization: roughness must be finite, positive and below 10 m.");
    if (!scene.getClimate())
        throw runtime_error("Voxelization: a readable Climate file is required for the wind profile.");

    log << "Voxelization: triangulating CitySim geometry..." << endl;
    SceneGeometry geometry=collectGeometry(scene,options);
    const size_t nx=max<size_t>(
                1,static_cast<size_t>(ceil(
                    (geometry.maximum.x-geometry.minimum.x)/options.deltaX)));
    const size_t ny=max<size_t>(
                1,static_cast<size_t>(ceil(
                    (geometry.maximum.y-geometry.minimum.y)/options.deltaY)));
    const size_t nz=max<size_t>(
                1,static_cast<size_t>(ceil(
                    (geometry.maximum.z-geometry.minimum.z)/options.deltaZ)));
    const size_t cells=checkedCellCount(nx,ny,nz,options.maximumCellCount);
    VoxelGrid grid(geometry.minimum.x,geometry.minimum.y,geometry.minimum.z,
                   options.deltaX,options.deltaY,options.deltaZ,nx,ny,nz);
    const int threadCount=effectiveThreadCount(options);
    const unsigned int windDay=scene.getBeginDay();
    const unsigned int windHour=1;
    grid.setLogLawWindProfile(
                scene.getClimate()->getWindSpeed(windDay,windHour),
                scene.getClimate()->getWindDirection(windDay,windHour),
                options.roughnessLength,10.);

    log << "Voxelization: AABB ["
        << geometry.minimum.x << "," << geometry.minimum.y << "," << geometry.minimum.z
        << "] -> [" << geometry.maximum.x << "," << geometry.maximum.y << "," << geometry.maximum.z
        << "], grid " << nx << "x" << ny << "x" << nz << " = " << cells
        << " cells, threads=" << threadCount << "." << endl;
    log << "Voxelization: neutral log-law wind at climate day=" << windDay
        << ", hour=" << windHour
        << ": U10=" << grid.referenceWindSpeed()
        << " m/s, direction=" << grid.referenceWindDirection()
        << " deg meteorological, roughness=" << grid.roughnessLength()
        << " m." << endl;

    const GroundReconstructionStats groundStats=
            voxelizeGround(geometry,grid,options,threadCount);
    log << "Voxelization: terrain seeds: explicit ground="
        << groundStats.explicitGroundColumns
        << " columns, building floors=" << groundStats.buildingFloorColumns
        << ", tree bases=" << groundStats.treeBaseColumns;
    if (groundStats.usedGlobalMinimum)
        log << "; no terrain evidence, used scene minimum z="
            << geometry.minimum.z;
    log << ". Missing columns were interpolated; the terrain covers the full domain."
        << endl;

    mergeBuildings(geometry,grid,options,threadCount);
    const bool hasBuildingShell=!geometry.buildings.empty() || !geometry.districtObstacles.empty();
    fillEnclosedBuildingCells(grid,hasBuildingShell,threadCount);
    voxelizeTrees(geometry.treeCanopies,grid,threadCount);
    log << "Voxelization: per-Tree wind attenuation assigned to "
        << grid.attenuatedTreeCellCount() << " TREE cells." << endl;

    const VoxelCounts counts=grid.counts();
    const double elapsed=chrono::duration<double>(
                chrono::steady_clock::now()-start).count();
    log << "Voxelization complete: AIR=" << counts.air
        << ", BUILDING=" << counts.building
        << ", GROUND=" << counts.ground
        << ", TREE=" << counts.tree
        << ", elapsed=" << elapsed << " s." << endl;
    return grid;
}
