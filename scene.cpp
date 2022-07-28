#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <limits>

#include "scene.h"
#include "district.h"
#include "plant.h"
#include "models.h"

#include "RENIndexedFaceSet.h" // for the attempt of triangulation

#ifndef DEBUG
  #include <omp.h> // for the OpenMP parallelisation
#endif

#ifdef FMI
// helper function for printing debug info
void importlogger(jm_callbacks* c, jm_string module, jm_log_level_enu_t log_level, jm_string message)
{
    printf("module = %s, log level = %s: %s\n", module, jm_log_level_to_string(log_level), message);
}
#endif

// *** Scene class, CitySim  *** //
// *** jerome.kaempf@epfl.ch *** //

Scene::Scene(string inputFile, string climateFile):inputFile(inputFile),logStream(std::cout.rdbuf()) {

    // loads the climate file
    if (!climateFile.empty()) {
        pClimate = new Climate(climateFile);

        // shows the information about the location
        logStream << "Location: " << pClimate->getLocation();
        logStream << "\t(Latitude: " << pClimate->getLatitudeN() << " N, Longitude: " << pClimate->getLongitudeE();
        logStream << " E, Altitude: " << pClimate->getAltitude() << " m, Meridian: " << pClimate->getMeridian() << " h)" << endl << flush;

        // gets the information from the climate file
        float latN = pClimate->getLatitudeN();
        float longE = pClimate->getLongitudeE();
        float merE = pClimate->getMeridian()*360.f/24.f;
        float northW = 0.f;

        // initilisation of the location for the sun and the scene
        SKYSiteLocation location(GENAngle::Degrees(latN), GENAngle::Degrees(longE), GENAngle::Degrees(merE), GENAngle::Degrees(northW));
        pSun = new SKYSun(location);
        scene.SetLocation(location);
    }

    // rapid initialisation of the parameters (to avoid segmentation faults)
    mNbReflections = 0;
    lv.assign(tregenzaSky.PatchCount()/2, 0.f);
    groundRadiance = 0.f;

}

void Scene::computeViewFactors() {

    // calculates the view factors of the scene
    // N.B.: the view factor calculation starts the direct, diffuse and daylight calculations in sequence
    // the direct calculation needs the Site Location in order to compute all sun positions
    v.CalculateViewFactors(scene);
    logStream << "View factors calculated." << endl << flush;

    // computes the projected solid angles of the vault (ground, sky and hemisphere) -> to get the SVF, GVF for each surface
    computeProjectedSolidAngles();

    // computes the sparse matrix for inter-reflections
    buildSparseMatrix();

}

void Scene::computeProjectedSolidAngles() {

    // an hemisphere is attached to each surface of the scene, this method
    // computes the projected solid angles of the vault (ground, sky and hemisphere) -> to get the SVF, GVF for each surface
    vector<double> sky(scene.SurfaceCount(), 0.), ground(scene.SurfaceCount(), 0.), hemisphere(scene.SurfaceCount(), 0.); // the seen projected hemisphere from and on the surface itself
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
        // diffuse part
        for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).SWViewFactors().GetVFs();
            factors!=scene.GetSurface(surfaceIndex).SWViewFactors().GetLastVF();
            ++factors) // loop on the patches that are non-zero
        {
            // unobstructed part of the Tregenza patch -> meaning not obstructed by a surface
            // the factors->unobstructed contains the patchSolidAngle and the cos of the angles between the surface normal and the center of patch
            if (factors->unobstructed > 0.f) {
                if ( factors->patchNo < static_cast<int>(tregenzaSky.PatchCount()/2) ) {
                    // diffuse part sky
                    sky[surfaceIndex] += factors->unobstructed;
                }
                else {
                    // diffuse part ground
                    ground[surfaceIndex] += factors->unobstructed;
                }
            }
            // the sum of the hemisphere (around Pi)
            hemisphere[surfaceIndex] += factors->unobstructed + factors->obstructed;
        }
        // saves all the projected solid angles
        ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->setProjectedSolidAngle(sky[surfaceIndex], ground[surfaceIndex], hemisphere[surfaceIndex]);
    }

    return;

}

void Scene::buildSparseMatrix() {

    logStream << "Building inter-reflection matrix." << endl << flush;
    logStream << "Number of surfaces: " << scene.SurfaceCount() << endl << flush;
    logStream << "Number of reflections: " << mNbReflections << endl << flush;

    // initialise the parameters
    Ai.clear();
    Aj.clear();
    An.clear();

    #ifdef DEBUG
    ostringstream saveVault;
    // for each patch
    for (unsigned int i=0; i<tregenzaSky.PatchCount(); ++i)
    {
        saveVault << i << "\t" << fmod(tregenzaSky.GetPatch(i)->centroid().Azimuth().degrees()+360.,360.) << "\t"
                  << tregenzaSky.GetPatch(i)->centroid().Altitude().degrees() << "\t"
                  << tregenzaSky.GetPatch(i)->solidAngle() << endl << flush;

        // for each Cell in the patch
        for (unsigned int j=0; j<tregenzaSky.GetPatch(i)->cellCount(); ++j)
        {
                saveVault << i << "\t" << j << "\t" << tregenzaSky.GetPatch(i)->getCell(j).direction().Azimuth().degrees() << "\t"
                          << tregenzaSky.GetPatch(i)->getCell(j).direction().Altitude().degrees() << "\t"
                          << tregenzaSky.GetPatch(i)->getCell(j).solidAngle() << endl << flush;
        }
    }
    save(string("vault.txt"),saveVault);
    #endif

    // only interesting when considering reflections
    if ( mNbReflections == 0 ) return;

    // loop on all surfaces to get the values in the matrices
    for (unsigned int surfaceIndex = 0; surfaceIndex < scene.SurfaceCount(); surfaceIndex++) {   // selection of a surface

            // the surface index is the row number
//            logStream << "Surface index: " << surfaceIndex << endl << flush;
//            logStream << "Surface ID: " << ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getId() << endl << flush;

            // add the reference in the Ai vector for surfaceIndex (corresponding line i)
            Ai.push_back( An.size() );

            // stores the end of the actual vector of row indices
            unsigned int endOfRowIndex = Aj.size();

            // going through the Tregenza patches
            for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).SWViewFactors().GetVFs();
                factors!=scene.GetSurface(surfaceIndex).SWViewFactors().GetLastVF();
                ++factors)
            {

                // if the obstructed fraction is null, nothing is taken into account (many values equal to zero come out from SRA, when unobstructed > 0.f)
                if ( !(factors->obstructed > 0.f) ) continue;

                // is the main obstructing surface already in the set?
                vector<unsigned int>::iterator it;
                it = find(Aj.begin()+endOfRowIndex, Aj.end(), factors->mainobstructing);
                if ( it == Aj.end() ) { // not found in the set
                    An.push_back( ( ((Surface*)(scene.GetSurface(factors->mainobstructing).SurfaceDelegate()))->getShortWaveReflectance() / M_PI )*(factors->obstructed) );
                    Aj.push_back( factors->mainobstructing );
                }
                else {
                    An[distance(Aj.begin(),it)] += ( ((Surface*)(scene.GetSurface(factors->mainobstructing).SurfaceDelegate()))->getShortWaveReflectance() / M_PI )*(factors->obstructed);
                }

//                logStream << "mainobstructing: " << factors->mainobstructing << "\tobstructed: " << factors->obstructed << "\tfor patch: " << factors->patchNo << endl << flush;
//                logStream << "An: " << An << endl << flush;
//                logStream << "Aj: " << Aj << endl << flush;
//                logStream << "Ai: " << Ai << endl << flush;

            }
    }
    // shows the matrix in CRS format
    if (An.empty()) {
        logStream << "Empty Compressed Row Storage inter-reflection Sparse Matrix." << endl << flush;
    }
    else {
        //logStream << "An: " << An << endl << flush;
        //logStream << "Aj: " << Aj << endl << flush;
        //logStream << "Ai: " << Ai << endl << flush;
        logStream << "Compressed Row Storage inter-reflections Sparse Matrix ready." << endl << flush;
        logStream << "Surfaces: " << getnAi() << "\tNon-zero elements: " << An.size() << "\tAverage non-zero elements / surface: " << static_cast<float>(An.size())/getnAi() << endl << flush;
    }

    return;

}

void Scene::showViewFactors() {

    // output on the screen of the view factors

    for (DATASurfaceIterator it=scene.GetAllSurfaces();
        !it.isAtEnd();
        ++it)
    {
        // selection of a surface it

        for (DATAViewFactorSetSparse::const_iterator factors=it->SWViewFactors().GetVFs();
            factors!=it->SWViewFactors().GetLastVF();
            ++factors)
        {
            // for that surface it, output of the informations in the sparse format

            logStream << "Unobstructed factor to patch " << factors->patchNo << ": " << factors->unobstructed << "\n";
            logStream << "Obstructed factor to the same patch: " << factors->obstructed << endl << flush;
            logStream << "Main obstructing surface: " << factors->mainobstructing << endl << flush;

        }
        pSun->SetDay(30); //absolute day from the start of the year [1,365]
        pSun->SetClockTime(10.5); // clock time [0,24[
        logStream << "Insolation factor: " << it->InsolationFactors().GetInsolationFactor(*pSun) << "\n" << flush;
    }

}

void Scene::exportRadFile(string radFile,bool triangulated) {

    // keeps the name of the inputFile if radFile empty
    if (radFile.empty()) radFile=inputFile;

    // open the output file
    ofstream outputRad((radFile.substr(0,radFile.size()-4) + ".rad").c_str(), ios::binary);

    outputRad.setf(ios::fixed); // set fixed floating format
    outputRad.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    outputRad.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // test d'ouverture
    if (!outputRad.is_open()) throw string("Error opening file: " + (radFile.substr(0,radFile.size()-4) + ".rad"));

    // writing the material
    outputRad << "void plastic surface\n" << "0\n" << "0\n" << "5 .2 .2 .2 0 0\n" << endl;

    if (triangulated) {

        // get Indexed Triangulated Surfaces
        std::shared_ptr<RENIndexedFaceSet> pIndexedFaceSet = scene.IndexedFaceSet();
        const std::vector<unsigned int>& pointIndices=pIndexedFaceSet->GetIndices();

        // loop on triangles
        for (unsigned int i=0; i<pointIndices.size()/3; i++) {

            outputRad << "surface polygon SURFACE#" << i << endl;
            outputRad << "0\n" << "0\n";
            // 9 coordinates given
            outputRad << "9\n";
            for (unsigned int k=0; k<3; k++) outputRad << GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3]))[k] << " ";
            outputRad << endl;
            for (unsigned int k=0; k<3; k++) outputRad << GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3+1]))[k] << " ";
            outputRad << endl;
            for (unsigned int k=0; k<3; k++) outputRad << GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3+2]))[k] << " ";
            outputRad << "\n" << endl;

        }
    }
    else {

        // loop on all surfaces
        for (DATASurfaceIterator it=scene.GetAllSurfaces();
            !it.isAtEnd();
            ++it)
        {

            outputRad << "surface polygon SURFACE#" << scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end()) << endl;
            outputRad << "0\n" << "0\n";

            // write the number of vertices for surface it
            outputRad << 3*(it->SurfaceDelegate()->vertexCount()) << endl;

            // get the surface vertices
            struct Vertices : public DATASurfaceDelegateABC::VertexVisitor
            {
                virtual void operator()(const GENPoint& v)
                {
                    vertices.push_back(v);
                }
                virtual ~Vertices() {};

                std::vector<GENPoint> vertices;
            } polygon;
            it->SurfaceDelegate()->sendVertices(polygon);

            // loop on the vertices to write them out
            for (unsigned int i=0; i<it->SurfaceDelegate()->vertexCount(); i++) {

                outputRad << "\t" << polygon.vertices[i][0] << " "
                                  << polygon.vertices[i][1] << " "
                                  << polygon.vertices[i][2] << endl;

            }
            outputRad << endl;

        }

    }

    outputRad.close();

}

void Scene::exportInpFile(string radFile, bool buildingsOnly) {

    // keeps the name of the inputFile if radFile empty
    if (radFile.empty()) radFile=inputFile;

    // open the output file
    ofstream outputInp((radFile.substr(0,radFile.size()-4) + ".inp").c_str(), ios::binary);

    outputInp.setf(ios::fixed); // set fixed floating format
    outputInp.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    outputInp.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // test d'ouverture
    if (!outputInp.is_open()) throw string("Error opening file: " + (radFile.substr(0,radFile.size()-4) + ".inp"));

    if (buildingsOnly) {
        // loop on all surfaces
        for (DATASurfaceIterator it=scene.GetBuildingSurfaces();
            !it.isAtEnd();
            ++it)
        {
            double gap = 0.05;

            // selection of a surface it
            //outputInp << "#Surface number: " << scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end()) << endl;

            GENPoint centroid = it->Centroid() + (it->Normal())*gap;

            outputInp << centroid[0] << "\t" << centroid[1] << "\t" << centroid[2] << "\t";
            outputInp << it->Normal()[0] << "\t" << it->Normal()[1] << "\t" << it->Normal()[2] << endl;
        }
    }
    else {
        // loop on all surfaces
        for (DATASurfaceIterator it=scene.GetAllSurfaces();
            !it.isAtEnd();
            ++it)
        {
            double gap = 0.05;

            // selection of a surface it
            //outputInp << "#Surface number: " << scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end()) << endl;

            GENPoint centroid = it->Centroid() + (it->Normal())*gap;

            outputInp << centroid[0] << "\t" << centroid[1] << "\t" << centroid[2] << "\t";
            outputInp << it->Normal()[0] << "\t" << it->Normal()[1] << "\t" << it->Normal()[2] << endl;
        }
    }

    outputInp.close();

}

void Scene::exportDXF(string fileName) {

    if(fileName=="") fileName=inputFile.substr(0,inputFile.size()-4) + ".dxf";
    else if(fileName.substr(fileName.size()-4,4)!=".dxf") fileName.append(".dxf");

    // open the output file
    ofstream outputDxf(fileName.c_str());
    if (!outputDxf.good()) throw(string("Error creating model: " + fileName));

    // setting good precision for the points
    outputDxf.setf(ios::fixed); // set fixed floating format
    outputDxf.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    outputDxf.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // writing the HEADER
    outputDxf << "999\n" << "CitySim DXF" << endl;
    outputDxf << "0\n" << "SECTION" << endl;
    outputDxf << "2\n" << "HEADER" << endl;
    outputDxf << "9\n" << "$ACADVER" << endl;
    outputDxf << "1\n" << "AC1006" << endl;
    outputDxf << "9\n" << "$INSUNITS" << endl;
    outputDxf << "70\n" << "6" << endl;
    outputDxf << "0\n" << "ENDSEC" << endl;

    // section of the ENTITIES
    outputDxf << " 0\n" << "SECTION\n" << " 2\n" << "ENTITIES" << endl;

    // loop on the surfaces for the export
    for (DATASurfaceIterator it=scene.GetAllSurfaces();
        !it.isAtEnd();
        ++it)
    {

        outputDxf << " 0\n" << "3DFACE" << endl;

        outputDxf << "8" << endl; // now comes the layer, 1 is for wall and 0 is for the ground
        if (it->IsBuildingSurface()) { outputDxf << "1" << endl; }
        else { outputDxf << "0" << endl; }
        outputDxf << "62" << endl; // 62 stands for the color of the surface
        if (it->IsBuildingSurface()) { outputDxf << "1" << endl; } // red for the buildings
        else { outputDxf << "2" << endl; } // yellow for the ground

        // now the vertices...
        struct Vertices : public DATASurfaceDelegateABC::VertexVisitor
        {
            virtual void operator()(const GENPoint& v)
            {
                vertices.push_back(v);
            }
            virtual ~Vertices() {}

            std::vector<GENPoint> vertices;
        } polygon;
        it->SurfaceDelegate()->sendVertices(polygon);

        if ( polygon.vertices.size() > 4 ) {
            // proceed to the triangulation
            std::shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
            RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
            meshBuilder.AddSurface(polygon.vertices.begin(),polygon.vertices.end(),0);
            const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
            for (unsigned int i=0; i<indexedFaceSet->TriangleCount()*3; i+=3) {
                // write those three points and repeat the first one
                for (unsigned int j=0; j<4; ++j) {
                    outputDxf << " 1" << i << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+(j%3)]))[0] << "\n"
                              << " 2" << i << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+(j%3)]))[1] << "\n"
                              << " 3" << i << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+(j%3)]))[2] << endl;
                }
            }
        }
        else {
            // loop on the vertices to write them out
            for (unsigned int i=0; i<polygon.vertices.size(); i++) {

                outputDxf << " 1" << i << "\n" << polygon.vertices[i][0] << "\n"
                          << " 2" << i << "\n" << polygon.vertices[i][1] << "\n"
                          << " 3" << i << "\n" << polygon.vertices[i][2] << endl;
            }
            if ( polygon.vertices.size() == 3 ) { // to close the loop the the 3 vertices

                outputDxf << " 13\n" << polygon.vertices[0][0] << "\n"
                          << " 23\n" << polygon.vertices[0][1] << "\n"
                          << " 33\n" << polygon.vertices[0][2] << endl;
            }
        }
    }

    outputDxf << " 0\n" << "ENDSEC\n" << " 0\n" << "EOF" << endl;
    outputDxf.close();

}

void Scene::computeRadiance(const unsigned int& day, const float& Idh, const float& Ibn, const float& albedo) {

    // erases and prepares the new vector
    lv.assign( tregenzaSky.PatchCount()/2, 0.f);
    groundRadiance = 0.f;

    // computes the 145 Tregenza patches' Radiance
    if (sky.SetSkyConditions(Idh,Ibn,pSun->GetPosition().Altitude().radians(), pSun->GetPosition().Azimuth().radians(),day)) {
        // if true Idh is positive and the diffuse sky can be created
        double intHemisphere = 0.;
        for (unsigned int i=0; i<tregenzaSky.PatchCount()/2; i++) {
            lv[i] = sky.GetRelativeLuminance(tregenzaSky.getPatchCenterAltitude(i),tregenzaSky.getPatchCenterAzimuth(i));
            intHemisphere += lv[i]*tregenzaSky.GetPatch(i)->solidAngle()*sin(tregenzaSky.getPatchCenterAltitude(i));
        }
        if (intHemisphere > 0.) for (unsigned int i=0; i<tregenzaSky.PatchCount()/2; i++) lv[i]*=Idh/intHemisphere;

        // computes the ground radiance value
        // as a lambertian reflection of the total sky and sun irradiances on the horizontal plane
        for (unsigned int p = 0; p < tregenzaSky.PatchCount()/2; p++){
            groundRadiance += lv[p] * tregenzaSky.GetPatch(p)->formFactor( GENPoint::Cartesian(0.f,0.f,1.f) ); // projected solidAngle in the horizontal plane
        }
        groundRadiance += (Ibn * sin(pSun->GetPosition().Altitude().radians())); //add solar irradiance, transform to beam horizontal
        groundRadiance *= albedo/M_PI;
    }

    return;
}

void Scene::computeCumulativeRadiance(unsigned int beginDay, unsigned int endDay, float albedo) {

    // vector to save the cumulative radiance
    vector<float> cumLv(tregenzaSky.PatchCount()/2, 0.f);
    float cumGroundRadiance = 0.f;

    // debut de la simulation sur la periode consideree
    for (unsigned int day = beginDay; day<=endDay; ++day) {
        for (unsigned int hour = 1; hour <= 24; ++hour) {

            // initialisation of the sun (for the VFC)
            pSun->SetDay(day);
            if (pSun->SetClockTime1(hour)) { // if sun is up

                // computes the patches radiance without obstructions
                computeRadiance(day, pClimate->getIdh(day,hour), pClimate->getIbn(day,hour), albedo);

                // accumulate it in the cumulative vector and value
                for (size_t i = 0; i < tregenzaSky.PatchCount()/2; ++i) {
                    cumLv[i] += lv[i];
                }
                cumGroundRadiance += groundRadiance;

            }
        }
    }

    // save the result in lv and groundRadiance
    for (size_t i = 0; i < tregenzaSky.PatchCount()/2; ++i) { lv[i] = cumLv[i]; }
    groundRadiance = cumGroundRadiance;

    // save the data
    //save("cumulativeSkyRadiance.txt", lv);
    //save("cumulativeGroundRadiance.txt", groundRadiance);

    // save as well the vault data
//  ostringstream saveVault;
    // for each patch
//	for (size_t i=0; i<tregenzaSky.PatchCount(); ++i)
//	{
//	    saveVault << i << "\t" << fmod(tregenzaSky.GetPatch(i)->centroid().Azimuth().degrees()+360.f,360.f) << "\t"
//                  << tregenzaSky.GetPatch(i)->centroid().Altitude().degrees() << "\t"
//                  << tregenzaSky.GetPatch(i)->solidAngle() << endl;
//	}
//    save("skyVault.txt",saveVault);

    return;

}

void Scene::exportSkyAndGround(string radFile, float luminousEfficacy) {

    // create the rad file with reference to .cal file
    ofstream radOut(radFile.c_str(), ios::binary);
    if (!radOut.is_open()) throw(string("Error creating sky file: " + radFile));

    radOut << "#Sky and Ground file generated by CitySim (jerome.kaempf@kaemco.ch)"
           << "\n" << endl;

    // preparation of the sky
    radOut << "void brightfunc skyfunc" << endl
           << "2 skybright " << radFile.substr(0,radFile.size()-4) << ".cal" << endl
           << "0" << endl << "0\n" << endl
           << "skyfunc glow sky_glow" << endl
           << "0" << endl << "0" << endl << "4 1 1 1 0\n" << endl
           << "sky_glow source sky" << endl
           << "0" << endl << "0" << endl << "4 0 0 1 180\n" << endl
           << "void glow ground_glow\n"
           << "0\n0\n4 " << groundRadiance*luminousEfficacy << " " << groundRadiance*luminousEfficacy << " " << groundRadiance*luminousEfficacy << " 0\n" << endl
           << "ground_glow source ground"
           << "\n0\n0\n4 0 0 -1 180";

    radOut.close();

    // create the .cal file with the radiance values
    ofstream skyCalOut((radFile.substr(0,radFile.size()-4) + ".cal").c_str(), ios::binary);
    if (!skyCalOut.is_open()) throw(string("Error creating sky file: " + radFile.substr(0,radFile.size()-4) + ".cal"));

    skyCalOut << "skybright=";
    for (unsigned int j=0; j<tregenzaSky.getBands()-1; j++)
    {
        skyCalOut << "row" << j << "+";
    }
    skyCalOut << "row" << tregenzaSky.getBands()-1 << ";" << endl << endl;

    unsigned int counter = 0;
    for (unsigned int j=0; j<tregenzaSky.getBands()-1; j++)
    {
        // note first patch split into two parts - first part (> 0 deg) and last patch (<360)
        skyCalOut << "row" << j << "=if(and(alt-" << j*tregenzaSky.getDeltaAltitude()*180./M_PI << ", " << (j+1)*tregenzaSky.getDeltaAltitude()*180./M_PI << "-alt),";
        skyCalOut << "select(floor(0.5+az/" << tregenzaSky.getDeltaAzimuth(j)*180./M_PI << ")+1," << endl;

        for (unsigned int i=counter; i< counter + tregenzaSky.getPatchesPerBand(j); i++)
        {
            skyCalOut << "\t" << lv[i]*luminousEfficacy << "," << endl;
        }
        // rewrite the first one.
        skyCalOut << "\t" << lv[counter]*luminousEfficacy << "),0);" << endl << endl;
        counter += tregenzaSky.getPatchesPerBand(j);
    }

    // top patch.
    skyCalOut << "row" << tregenzaSky.getBands()-1 << "=if(alt-"<< 90.-(tregenzaSky.getDeltaAltitude()*180./M_PI / 2.)<< "," << lv[counter]*luminousEfficacy << ",0);"<< endl << endl;

    skyCalOut << "alt=asin(Dz)*180/PI;" << endl << endl;
    skyCalOut << "az=if(azi,azi,azi+360);" << endl;
    skyCalOut << "azi=atan2(Dx,Dy)*180/PI;" << endl << endl;

    skyCalOut.close();

    return;
}

void Scene::exportCumulativeRadiance()
{
    // computes the sky and ground file
    computeCumulativeRadiance();
    exportSkyAndGround((inputFile.substr(0,inputFile.size()-4) + "_cumSkyGrnd.rad"));

    // preparation of the cumulative suns file
    ofstream radOut((inputFile.substr(0,inputFile.size()-4) + "_cumSuns.rad").c_str(), ios::binary);
    if (!radOut.is_open()) throw string("Error creating the cumulative Suns file.");

    double sunLuminance = 0.;
    GENPoint sunPosition;
    for (unsigned int day = 1; day<=365; ++day) {
        for (unsigned int hour = 1; hour <= 24; ++hour) {

            // initialisation of the sun (for the VFC)
            pSun->SetDay(day);
            pSun->SetClockTime1(hour);
            if (pSun->SunUp()) { // if sun is up

                sunLuminance = pClimate->getIbn(day,hour) / pSun->getSolidAngle();

                radOut << "void light solar_" << day << "_" << hour << "\n"
                       << "0\n"
                       << "0\n"
                       << "3 " << sunLuminance << " " << sunLuminance << " " << sunLuminance << endl;

                sunPosition = pSun->GetPosition();

                radOut << "solar_" << day << "_" << hour << " source sun_" << day << "_" << hour << "\n"
                       << "0\n"
                       << "0\n"
                       << "4 " << sunPosition[0] << " "
                               << sunPosition[1] << " "
                               << sunPosition[2] << " "
                               << pSun->getAperture() << "\n"
                       << endl;
            }
        }
    }

}

// THE RADIANCE SCENE

Radscene::Radscene(string inputFile, string climateFile):Scene(inputFile, climateFile) {

    // reads the Radiance file description and put the facades in the surfaceVector
    logStream << "Reading Radiance file..." << endl << flush;

    string tampon;
    char tampon1[200];

    // chargement des donnees a afficher
    ifstream input1 (inputFile.c_str(), ios::binary);

    // test d'ouverture
    if (!input1.is_open()) throw string("Error opening Radiance file");


    // initialise le compteur de surfaces
    unsigned int count = 0;
    while (!input1.eof()) {

        input1 >> tampon; // first keyword

        if ( tampon[0] == '#' ) input1.getline(tampon1, 200,'\n');
        else if ( tampon == "void" ) { // material definition

            input1 >> tampon; // material type
            input1 >> tampon; // label of material
            input1 >> tampon; // 0
            input1 >> tampon; // 0
            input1 >> tampon; // number of parameters
            unsigned int number = atoi(tampon.c_str());

            for (unsigned int i=0; i<number; i++) input1 >> tampon;

        }
        else { // whatever the name of the polygon

            input1 >> tampon; // polygon
            if ( tampon == "polygon" ) { // it has to be a polygon

                if (input1.eof()) break;

                // create the polygon
                //surfaceVector.push_back(GENHandle<MySurfaceDelegate>(new MySurfaceDelegate()));
                GENHandle<Surface> surface(new Surface(count++,0.2f,0.f,0.f,0.f,0.f));

                input1 >> tampon; // label
                input1 >> tampon; // 0
                input1 >> tampon; // 0
                input1 >> tampon; // number of points (3x3)
                unsigned int number = atoi(tampon.c_str())/3;

                double x,y,z;
                for (unsigned int i=0; i<number; i++) {

                    input1 >> tampon;
                    x = atof(tampon.c_str());
                    input1 >> tampon;
                    y = atof(tampon.c_str());
                    input1 >> tampon;
                    z = atof(tampon.c_str());

                    logStream << "Point: " << x << ", " << y << ", " << z << endl << flush;

                    //surfaceVector.back()->pushVertex(x,y,z);
                    surface->pushVertex(x,y,z);

                }
                // computation of the surface normal (anti-clockwise -> positive)
                //surfaceVector.back()->computeNormal();
                surface->computeNormal();
                //logStream << "Normal computed: (" << surfaceVector.back()->normal()[0] << ",";
                //logStream << surfaceVector.back()->normal()[1] << "," << surfaceVector.back()->normal()[2] << ")" << endl << flush;
                logStream << "Normal computed: (" << surface->normal()[0] << ",";
                logStream << surface->normal()[1] << "," << surface->normal()[2] << ")" << endl << flush;
                // computation of the surface area
                //surfaceVector.back()->computeArea();
                surface->computeArea();
                //logStream << "Area computed: " << surfaceVector.back()->getArea() << endl << flush;
                logStream << "Area computed: " << surface->getArea() << endl << flush;
                if (surface->getArea() > 0) { // if non degenerate area
                    // adds the surface to the scene
                    scene.AddGroundSurface(surface); // no need of the daylighting model
                }
                else { logStream << "Null surface area not taken into account." << endl << flush; }
            }
            else {
                input1 >> tampon; // label of material
                input1 >> tampon; // number of ascii parameters
                unsigned int number = atoi(tampon.c_str());
                for (unsigned int i=0; i<number; i++) input1 >> tampon;
                input1 >> tampon; // 0
                input1 >> tampon; // number of parameters
                number = atoi(tampon.c_str());
                for (unsigned int i=0; i<number; i++) input1 >> tampon;
            }
        }
    }
    input1.close();

    // calculates the view factors of the scene
    // N.B.: the view factor calculation starts the direct, diffuse and daylight calculations in sequence
    // the direct calculation needs the Site Location in order to compute all sun positions
    computeViewFactors();

    // initialise the vectors for intermediary results
    irradiationDiffuseSky.assign(scene.SurfaceCount(),0.f);
    irradiationDiffuseGround.assign(scene.SurfaceCount(),0.f);
    irradiationBeam.assign(scene.SurfaceCount(),0.f);
    irradiationReflection.assign(scene.SurfaceCount(),0.f);

}

void Radscene::clearResults() {
    irradiationDiffuseSky.assign(scene.SurfaceCount(),0.f);
    irradiationDiffuseGround.assign(scene.SurfaceCount(),0.f);
    irradiationBeam.assign(scene.SurfaceCount(),0.f);
    irradiationReflection.assign(scene.SurfaceCount(),0.f);
    outData.str(""); // clearing the ostringstream
}

void Radscene::exportSWFile(string filename) {

  fstream output (filename.c_str(), ios::out | ios::binary);

  output.setf(ios::fixed); // set fixed floating format
  output.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
  output.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

  for (unsigned int i=0;i<irradiationBeam.size();++i) output << irradiationBeam[i]+irradiationDiffuseSky[i]+irradiationDiffuseGround[i]+irradiationReflection[i] << endl;

  output.close();

  return;

}

void Radscene::computeShortWave(unsigned int day, unsigned int hour) {

    // irradiationSW: irradiation without reflections, irradiationSWn: irradiation with reflections
    vector<float> irradiationSW(scene.SurfaceCount(), 0.f), irradiationSWn(irradiationSW);

    // initialisation of the sun (for the VFC)
    pSun->SetDay(day);
    if (pSun->SetClockTime1(hour)) { // if the day has started

        // gets Idh and Ibn from the climate file
        float Idh = pClimate->getIdh(day,hour);
        float Ibn = pClimate->getIbn(day,hour);

        // computes the patches radiance without and then with obstructions
        computeRadiance(day,Idh,Ibn);
        //save("cumRadiance_" + toString(day) + "_" + toString(hour) + ".txt", lv);

        // now the main loop on surfaces
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
            // diffuse part sky
            for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).SWViewFactors().GetVFs();
                factors!=scene.GetSurface(surfaceIndex).SWViewFactors().GetLastVF();
                ++factors) // loop on the patches that are non-zero
            {
                // the factors->unobstructed contains the patchSolidAngle and the cos of the angles between the surface normal and the center of patch
                if (factors->unobstructed > 0.f) {
                    if ( factors->patchNo < static_cast<int>(tregenzaSky.PatchCount()/2) ) {
                        // diffuse part sky
                        irradiationSW[surfaceIndex] += factors->unobstructed * lv[factors->patchNo];
                        irradiationDiffuseSky[surfaceIndex] += factors->unobstructed * lv[factors->patchNo];
                    }
                }
            }
            // diffuse part ground
            irradiationSW[surfaceIndex] += ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_ground() * groundRadiance;
            irradiationDiffuseGround[surfaceIndex] += ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_ground() * groundRadiance;
            // direct part
            double cosTheta = GEN::dot_product(scene.GetSurface(surfaceIndex).Normal(),pSun->GetPosition());
            if (cosTheta > 0.) {
                irradiationSW[surfaceIndex] += scene.GetSurface(surfaceIndex).InsolationFactors().GetInsolationFactor(*pSun) * Ibn * cosTheta; // fraction of the surface that is light by the sun
                irradiationBeam[surfaceIndex] += scene.GetSurface(surfaceIndex).InsolationFactors().GetInsolationFactor(*pSun) * Ibn * cosTheta; // fraction of the surface that is light by the sun
            }
        } // end the loop on the surfaces

        // adds the reflections on the surfaces
        irradiationSWn = irradiationSW;
        for (unsigned int r = 0; r < mNbReflections; ++r) {

            // using the CRS sparse matrix format, multiplication of the matrix by the former vector
            vector<float> irradiationSWnew(irradiationSW);
            #pragma omp parallel for schedule(dynamic)
            for (unsigned int i=0; i < getnAi(); ++i) { // loop on the number of elements
                for (unsigned int index=getAi(i); index < getAi(i+1); ++index) {
                    irradiationSWnew[i] += getAn(index) * irradiationSWn[getAj(index)];
                }
            }
            // saves the new irradiation vector in the current
            irradiationSWn = irradiationSWnew;

        }

    } // end if sun is up

    // computes what is due to inter-reflections
    for (unsigned int i=0; i<irradiationSWn.size(); ++i) irradiationReflection[i] = irradiationSWn[i] - irradiationSW[i];

    // end of the calculations, save in all surfaces the irradiation
    for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
    {
        // saves all the irradiation
        ((Surface*)(it->SurfaceDelegate()))->setShortWaveIrradiance(irradiationSWn[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

    }

    return;

}

void Radscene::compareWithRadianceExternalIrradiance(unsigned int day, unsigned int hour) {

    // test if the sun is up
    pSun->SetDay(day);
    if (pSun->SetClockTime1(hour)) {

        // loop on the surfaces to create the _mesh.inp
        ofstream outputInp((inputFile.substr(0,inputFile.size()-4) + "_mesh.inp").c_str(), ios::binary);
        if (!outputInp.is_open()) throw string("Error opening file: " + (inputFile.substr(0,inputFile.size()-4) + "_mesh.inp"));

        outputInp.setf(ios::fixed); // set fixed floating format
        outputInp.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
        outputInp.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

        // save the number of points per surface
        vector<unsigned int> numberOfPointsPerSurface;

        // loop on all surfaces
        for (DATASurfaceIterator it=scene.GetAllSurfaces();
            !it.isAtEnd();
            ++it)
        {

            double gap = 0.05;
            double maxDetectorArea = 1.0;

            // get the surface vertices
            struct Vertices : public DATASurfaceDelegateABC::VertexVisitor
            {
                virtual void operator()(const GENPoint& v)
                {
                    vertices.push_back(v);
                }
                virtual ~Vertices() {}

                std::vector<GENPoint> vertices;
            } polygon;
            it->SurfaceDelegate()->sendVertices(polygon);

            unsigned int numberOfPoints = 0;

            for (unsigned int i=0;i<it->SurfaceDelegate()->vertexCount();i++) {

                vector<GENPoint> triangle,grid;
                triangle.push_back(it->Centroid());
                triangle.push_back(polygon.vertices[i]);
                triangle.push_back(polygon.vertices[(i+1)%(it->SurfaceDelegate()->vertexCount())]);

                // creation of the grid
                gridTriangle(triangle,maxDetectorArea,grid);

                // save the number of vertices
                numberOfPoints += grid.size();

                // loop on the points in the grid
                for (unsigned int j=0;j<grid.size();j++) {

                    // output of the couple grid and normal
                    outputInp << grid[j][0]+it->Normal()[0]*gap << " " << grid[j][1]+it->Normal()[1]*gap << " " << grid[j][2]+it->Normal()[2]*gap << " ";
                    outputInp << it->Normal()[0] << " " << it->Normal()[1] << " " << it->Normal()[2] << endl;

                }
            }

            // put the number of points in the vector
            numberOfPointsPerSurface.push_back(numberOfPoints);
        }
        outputInp.close();

        // vectors to save the results
        vector<float> sunIrrad, skyIrrad, groundIrrad, reflectionIrrad;

        // export the sun, the sky and the ground in Radiance format
        exportSunRadFile(inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sun.rad",day,hour);
        exportSkyRadFile(inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sky.rad",day,hour);
        exportGroundRadFile(inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_ground.rad",day,hour);
        // export the .inp files for the direct and diffuse
        exportInpFile();

        int error = 0; // error handling for the system command

        // computes the direct component, start Radiance and read the results
        ostringstream command;
        command << "oconv " << inputFile << " " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day)
                << "_" << toString(hour) << "_sun.rad > " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_"
                << toString(hour) << "_sun.oct";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start oconv (.oct)"));
        command.str("");
        command << "rtrace -w -I -h -dj 0 -ds 0 -dt 0 -dc 0.6 -dr 0 -dp 0 -ab 0 -aa 0 "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << "_sun.oct"
                << " < " << inputFile.substr(0,inputFile.size()-4) << "_mesh.inp | rcalc -e '$1=0.265*$1+0.67*$2+0.065*$3' > "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << "_sun.out";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start rtrace (.out)"));

        // computes the sky component
        command.str("");
        command << "oconv " << inputFile << " " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day)
                << "_" << toString(hour) << "_sky.rad > " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_"
                << toString(hour) << "_sky.oct";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start oconv (.oct)"));
        command.str("");
        command << "rtrace -w -I -h -dj 0 -ds 0 -dt 0 -dc 0.6 -dr 0 -dp 0 -ab 1 -ad 2048 -as 512 -aa 0 "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << "_sky.oct"
                << " < " << inputFile.substr(0,inputFile.size()-4) << ".inp | rcalc -e '$1=0.265*$1+0.67*$2+0.065*$3' > "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << "_sky.out";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start rtrace (.out)"));

        // computes the ground component
        command.str("");
        command << "oconv " << inputFile << " " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day)
                << "_" << toString(hour) << "_ground.rad > " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_"
                << toString(hour) << "_ground.oct";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start oconv (.oct)"));
        command.str("");
        command << "rtrace -w -I -h -dj 0 -ds 0 -dt 0 -dc 0.6 -dr 0 -dp 0 -ab 1 -ad 2048 -as 512 -aa 0 "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << "_ground.oct"
                << " < " << inputFile.substr(0,inputFile.size()-4) << ".inp  | rcalc -e '$1=0.265*$1+0.67*$2+0.065*$3' > "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << "_ground.out";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start rtrace (.out)"));

        // computes the inter-reflection
        command.str("");
        command << "oconv " << inputFile << " " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day)
                << "_" << toString(hour) << "_sun.rad" << " " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day)
                << "_" << toString(hour) << "_sky.rad" << " " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day)
                << "_" << toString(hour) << "_ground.rad > " << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_"
                << toString(hour) << ".oct";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start oconv (.oct)"));
        command.str("");
        command << "rtrace -w -I -h -dj 0 -ds 0 -dt 0 -dc 0.6 -dr 0 -dp 0 -ab " << mNbReflections+1 << " -ad 2048 -as 512 -aa 0 "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << ".oct"
                << " < " << inputFile.substr(0,inputFile.size()-4) << "_mesh.inp  | rcalc -e '$1=0.265*$1+0.67*$2+0.065*$3' > "
                << inputFile.substr(0,inputFile.size()-4) << "_" << toString(day) << "_" << toString(hour) << "_reflection.out";
        error = system(command.str().c_str());
        if (error!=0) throw(string("Cannot start rtrace (.out)"));

        // lecture des fichiers de sorties
        readResults(inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sun.out",sunIrrad); // sun
        readResults(inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sky.out",skyIrrad); // sky
        readResults(inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_ground.out",groundIrrad); // ground
        readResults(inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_reflection.out",reflectionIrrad); // sky

        // average results on each surface
        averageIrradianceOnSurfaces(numberOfPointsPerSurface,sunIrrad);
        averageIrradianceOnSurfaces(numberOfPointsPerSurface,reflectionIrrad);

        for (unsigned int i=0;i<scene.SurfaceCount();++i) {
            outData << day << "\t" << hour << "\t" << i << "\t"
                    << irradiationBeam[i]+irradiationDiffuseSky[i]+irradiationDiffuseGround[i]+irradiationReflection[i] << "\t"
                    << irradiationBeam[i] << "\t"
                    << irradiationDiffuseSky[i] << "\t"
                    << irradiationDiffuseGround[i] << "\t"
                    << reflectionIrrad[i] << "\t"
                    << sunIrrad[i] << "\t"
                    << skyIrrad[i] << "\t"
                    << groundIrrad[i] << "\t"
                    << scene.GetSurface(i).Area() << "\n";
        }

        // effacement des fichiers
        //remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sun.rad").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sun.oct").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sun.out").c_str());
        //remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sky.rad").c_str());
        //remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sky.cal").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sky.oct").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_sky.out").c_str());
        //remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_ground.rad").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_ground.oct").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_ground.out").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + ".oct").c_str());
        remove((inputFile.substr(0,inputFile.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_reflection.out").c_str());
    }
}

void Radscene::withRadianceClimate(string climateFile) {

    // exports the RAD file with which the calculation is made
    exportRadFile(inputFile.substr(0,inputFile.size()-4) + "_converted.rad");

    // loop on the surfaces to create the _mesh.inp
    ofstream outputInp((inputFile.substr(0,inputFile.size()-4) + "_mesh.inp").c_str(), ios::binary);
    if (!outputInp.is_open()) throw string("Error opening file: " + (inputFile.substr(0,inputFile.size()-4) + "_mesh.inp"));

    outputInp.setf(ios::fixed); // set fixed floating format
    outputInp.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    outputInp.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // save the number of points per surface
    vector<unsigned int> numberOfPointsPerSurface;

    // loop on all surfaces
    for (DATASurfaceIterator it=scene.GetAllSurfaces();
        !it.isAtEnd();
        ++it)
    {

        double gap = 0.05;
        double maxDetectorArea = 1.0;

        // get the surface vertices
        struct Vertices : public DATASurfaceDelegateABC::VertexVisitor
        {
            virtual void operator()(const GENPoint& v)
            {
                vertices.push_back(v);
            }
            virtual ~Vertices() {};

            std::vector<GENPoint> vertices;
        } polygon;
        it->SurfaceDelegate()->sendVertices(polygon);

        unsigned int numberOfPoints = 0;

        if (polygon.vertices.size() == 3) { // triangle

            vector<GENPoint> triangle,grid;
            triangle.push_back(polygon.vertices[0]);
            triangle.push_back(polygon.vertices[1]);
            triangle.push_back(polygon.vertices[2]);

            // creation of the grid
            gridTriangle(triangle,maxDetectorArea,grid);

            // save the number of vertices
            numberOfPoints += grid.size();

            // loop on the points in the grid
            for (unsigned int j=0;j<grid.size();j++) {

                // output of the couple grid and normal
                outputInp << grid[j][0]+it->Normal()[0]*gap << " " << grid[j][1]+it->Normal()[1]*gap << " " << grid[j][2]+it->Normal()[2]*gap << " ";
                outputInp << it->Normal()[0] << " " << it->Normal()[1] << " " << it->Normal()[2] << endl;

            }
        }
        else { // a triangulation would be nice
            throw(string("Error not a triangle!!"));
            for (unsigned int i=0;i<it->SurfaceDelegate()->vertexCount();i++) {

                vector<GENPoint> triangle,grid;
                triangle.push_back(it->Centroid());
                triangle.push_back(polygon.vertices[i]);
                triangle.push_back(polygon.vertices[(i+1)%(it->SurfaceDelegate()->vertexCount())]);

                // creation of the grid
                gridTriangle(triangle,maxDetectorArea,grid);

                // save the number of vertices
                numberOfPoints += grid.size();

                // loop on the points in the grid
                for (unsigned int j=0;j<grid.size();j++) {

                    // output of the couple grid and normal
                    outputInp << grid[j][0]+it->Normal()[0]*gap << " " << grid[j][1]+it->Normal()[1]*gap << " " << grid[j][2]+it->Normal()[2]*gap << " ";
                    outputInp << it->Normal()[0] << " " << it->Normal()[1] << " " << it->Normal()[2] << endl;

                }
            }
        }
        // put the number of points in the vector
        numberOfPointsPerSurface.push_back(numberOfPoints);
    }
    outputInp.close();

    int error = 0; // error handling for the system command

    // start Radiance, which directly writes the results in the file
    ostringstream command;
    command << "oconv " << inputFile.substr(0,inputFile.size()-4) << "_converted.rad " << climateFile << " > " << inputFile.substr(0,inputFile.size()-4) << ".oct";
    error = system(command.str().c_str());
    if (error!=0) throw(string("Cannot start oconv (.oct)"));
    command.str("");

    command << "rtrace -w -I -h -dj 0 -ds 0 -dt 0 -dc 0.6 -dr 0 -dp 0 -ab 2 -ad 2048 -as 512 -aa 0 "
            << inputFile.substr(0,inputFile.size()-4) << ".oct"
            << " < " << inputFile.substr(0,inputFile.size()-4) << "_mesh.inp  | rcalc -e '$1=0.265*$1+0.67*$2+0.065*$3' > "
            << inputFile.substr(0,inputFile.size()-4) << "_mesh.out";
    error = system(command.str().c_str());
    if (error!=0) throw(string("Cannot start rtrace (.out)"));

    // lecture des fichiers de sorties
    vector<float> irrad;
    readResults(inputFile.substr(0,inputFile.size()-4) + "_mesh.out",irrad); // sun

    // average results on each surface
    averageIrradianceOnSurfaces(numberOfPointsPerSurface,irrad);

    for (unsigned int i=0;i<scene.SurfaceCount();i++) {
        outData << irrad[i] << "\n";
    }
    outData << endl;

    // write the results in a .out file
    remove((inputFile.substr(0,inputFile.size()-4) + ".out").c_str());
    writeResults(inputFile.substr(0,inputFile.size()-4) + ".out" );

}

void Radscene::exportSunRadFile(string sunRadFile,unsigned int day,unsigned int hour) {

    // gets Ibn from the climate file
    float Ibn = pClimate->getIbn(day,hour);

    // initilisation of the sun (for the VFC)
    pSun->SetDay(day);
    pSun->SetClockTime1(hour);
    if (pSun->SunUp()) {

        // compute the solar radiance from the Ibn
        double solarRadiance= Ibn/pSun->getSolidAngle();

        // ouverture et ecriture du fichier
        ofstream sunFile(sunRadFile.c_str(), ios::binary);

        sunFile << "#Sun file generated by CitySim (jerome.kaempf@epfl.ch)\n"
                << "#Sun azimuth: " << pSun->GetPosition().Azimuth().degrees() << "\taltitude: " << pSun->GetPosition().Altitude().degrees()
                << "\tIbn: " << Ibn << "\n" << endl;

        sunFile << "void light solar" << endl;
        sunFile << "0" << endl << "0" << endl;
        sunFile << "3 " << solarRadiance << " " << solarRadiance << " " << solarRadiance << "\n" << endl;

        sunFile << "solar source sun" << endl;
        sunFile << "0" << endl << "0" << endl;
        sunFile << "4 " << pSun->GetPosition()[0] << " " << pSun->GetPosition()[1] << " " << pSun->GetPosition()[2] << " " << pSun->getAperture() << endl;

        sunFile.close();
    }
}

void Radscene::exportSkyRadFile(string skyRadFile,unsigned int day,unsigned int hour) {

    pSun->SetDay(day);
    pSun->SetClockTime1(hour);
    if (pSun->SunUp()) {

        // gets Idh and Ibn from the climate file
        float Idh = pClimate->getIdh(day,hour);
        float Ibn = pClimate->getIbn(day,hour);

        /// create the rad file with reference to .cal file

        ofstream skyRadOut(skyRadFile.c_str(), ios::binary);
        if (!skyRadOut.is_open()) throw(string("Error creating sky file: " + skyRadFile));

        skyRadOut << "#Sky file generated by CitySim (jerome.kaempf@epfl.ch)"
                  << "\n" << endl;

        skyRadOut << "void brightfunc skyfunc" << endl;
        skyRadOut << "2 skybright " << skyRadFile.substr(0,skyRadFile.size()-4) << ".cal" << endl;
        skyRadOut << "0" << endl << "0" << endl;
        skyRadOut << "skyfunc glow sky_glow" << endl;
        skyRadOut << "0" << endl << "0" << endl << "4 1 1 1 0" << endl;
        skyRadOut << "sky_glow source sky" << endl;
        skyRadOut << "0" << endl << "0" << endl << "4 0 0 1 180" << endl;

        skyRadOut.close();

        /// create the .cal file with the luminances values

        // get the 145 sky radiances and ground radiance
        computeRadiance(day,Idh,Ibn);

        // output the radiances
        ofstream skyCalOut((skyRadFile.substr(0,skyRadFile.size()-4) + ".cal").c_str(), ios::binary);
        if (!skyCalOut.is_open()) throw(string("Error creating sky file: " + skyRadFile.substr(0,skyRadFile.size()-4) + ".cal"));

        skyCalOut << "skybright=";
        for (unsigned int j=0; j<tregenzaSky.getBands()-1; j++)
        {
            skyCalOut << "row" << j << "+";
        }
        skyCalOut << "row" << tregenzaSky.getBands()-1 << ";" << endl << endl;

        unsigned int counter = 0;
        for (unsigned int j=0; j<tregenzaSky.getBands()-1; j++)
        {
            // note first patch split into two parts - first part (> 0 deg) and last patch (<360)
            skyCalOut << "row" << j << "=if(and(alt-" << j*tregenzaSky.getDeltaAltitude()*180./M_PI << ", " << (j+1)*tregenzaSky.getDeltaAltitude()*180./M_PI << "-alt),";
            skyCalOut << "select(floor(0.5+az/" << tregenzaSky.getDeltaAzimuth(j)*180./M_PI << ")+1," << endl;

            for (unsigned int i=counter; i< counter + tregenzaSky.getPatchesPerBand(j); i++)
            {
                skyCalOut << "\t" << lv[i] << "," << endl;
            }
            // rewrite the first one.
            skyCalOut << "\t" << lv[counter] << "),0);" << endl << endl;
            counter += tregenzaSky.getPatchesPerBand(j);
        }

        // top patch.
        skyCalOut << "row" << tregenzaSky.getBands()-1 << "=if(alt-"<< 90.-(tregenzaSky.getDeltaAltitude()*180./M_PI / 2.)<< "," << lv[counter] << ",0);"<< endl << endl;

        skyCalOut << "alt=asin(Dz)*180/PI;" << endl << endl;
        skyCalOut << "az=if(azi,azi,azi+360);" << endl;
        skyCalOut << "azi=atan2(Dx,Dy)*180/PI;" << endl << endl;

        skyCalOut.close();
    }
}

void Radscene::exportGroundRadFile(string groundRadFile,unsigned int day, unsigned int hour) {

    // gets Idh and Ibn from the climate file
    float Idh = pClimate->getIdh(day,hour);
    float Ibn = pClimate->getIbn(day,hour);

    computeRadiance(day,Idh,Ibn);

    // create the rad file with reference to .cal file
    ofstream groundRadOut(groundRadFile.c_str(), ios::binary);
    if (!groundRadOut.is_open()) throw(string("Error creating ground file: " + groundRadFile));

    groundRadOut << "void glow ground_glow\n"
                 << "0\n0\n4 " << groundRadiance << " " << groundRadiance << " " << groundRadiance << " 0\n"
                 << "\nground_glow source ground"
                 << "\n0\n0\n4 0 0 -1 180";
    groundRadOut.close();

}

void Radscene::readResults(string filename, vector<float> &results) {

    float mesure;
    string tampon;

    // chargement des donnees a afficher
    ifstream input1 (filename.c_str(), ios::binary | ios::in );

    // test d'ouverture

    if (!input1.is_open()) throw (string("Error opening: " + filename));

    results.clear();

    do {

        input1 >> tampon;

        if (tampon == "") break;

        sscanf(tampon.c_str(), "%f", &mesure);

        results.push_back(mesure);

        tampon = "";

    } while (!input1.eof());

    input1.close();

}

void Radscene::writeHeader(string filename) {

    // writes all outputs
    ofstream outputInFile(filename.c_str(), ios::binary | ios::app );

    // test d'ouverture
    if (!outputInFile.is_open()) throw(string("Error opening: " + filename));

    outputInFile << "#day\thour\tsurface\twithInter\tsun\tsky\tground\twithInter\tsun\tsky\tground\n";

    outputInFile.flush();
    outputInFile.close();

}

void Radscene::writeResults(string filename) {

    // writes all outputs
    ofstream outputInFile(filename.c_str(), ios::binary | ios::app );

    // test d'ouverture
    if (!outputInFile.is_open()) throw(string("Error opening: " + filename));

    outputInFile << outData.str();

    outputInFile.flush();
    outputInFile.close();

}

void Radscene::gridTriangle(vector<GENPoint> triangle, double maxDetectorArea, vector<GENPoint> &grid) {

    // decomposition du triangle en 4 sous-triangles et ainsi de suite

    double l1, l2, l3, area;

    vector<double> d1;
    vector<double> d2;

    vector<GENPoint> subdiv, subdiv2, triangle2;

    // calcul de l1, l2 et l3

    l1 = sqrt( pow( (triangle[1][0]-triangle[0][0]), 2.)
                + pow( (triangle[1][1]-triangle[0][1]), 2.)
                + pow( (triangle[1][2]-triangle[0][2]), 2.) );

    l2 = sqrt( pow( (triangle[2][0]-triangle[1][0]), 2.)
                + pow( (triangle[2][1]-triangle[1][1]), 2.)
                + pow( (triangle[2][2]-triangle[1][2]), 2.) );

    l3 = sqrt( pow( (triangle[0][0]-triangle[2][0]), 2.)
                + pow( (triangle[0][1]-triangle[2][1]), 2.)
                + pow( (triangle[0][2]-triangle[2][2]), 2.) );

//    logStream << "l1: " << l1 << "\tl2: " << l2 << "\tl3: " << l3 << endl << flush;

    // calcul de l'aire

    area = 0.25*sqrt( (l1 + l2 + l3)*(-l1 + l2 + l3)*(l1 - l2 + l3)*(l1 + l2 - l3) );

 //   logStream << "area: " << area << endl << flush;

    if (area < 1.e-6) return;

     // calcul de d1 et d2

    d1.push_back( (triangle[1][0]-triangle[0][0]) / l1 );
    d1.push_back( (triangle[1][1]-triangle[0][1]) / l1 );
    d1.push_back( (triangle[1][2]-triangle[0][2]) / l1 );

    d2.push_back( (triangle[2][0]-triangle[0][0]) / l3 );
    d2.push_back( (triangle[2][1]-triangle[0][1]) / l3 );
    d2.push_back( (triangle[2][2]-triangle[0][2]) / l3 );

    // distinction des cas

    if ( area <= maxDetectorArea ) {

            grid.clear();
            grid.push_back( centroid(triangle) );

    }
    else {

        unsigned int level;

        // preparation de subdiv: les 3 points du triangle

        subdiv.push_back(triangle[0]);
        subdiv.push_back(triangle[1]);
        subdiv.push_back(triangle[2]);

        // calcul du niveau de division a obtenir

        level = (unsigned int) int( ceil( log( area / maxDetectorArea )/log(4.0) ) );

//        logStream << "niveau: " << niveau << endl << flush;

        for (unsigned int i=0; i < level; i++) {

            grid.clear();
            subdiv2.clear();

            area /= 4.; // subdivision de la surface par point

            for (unsigned int j=0; j < (subdiv.size()/3); j++) {

                triangle2.clear();

                triangle2.push_back( middlePoint(subdiv[j*3],   subdiv[j*3+1]) );
                triangle2.push_back( middlePoint(subdiv[j*3+1], subdiv[j*3+2]) );
                triangle2.push_back( middlePoint(subdiv[j*3+2], subdiv[j*3]) );
                grid.push_back( centroid(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

                triangle2.clear();

                triangle2.push_back( subdiv[j*3] );
                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+1]) );
                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+2]) );
                grid.push_back( centroid(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

                triangle2.clear();

                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+1]) );
                triangle2.push_back( subdiv[j*3+1] );
                triangle2.push_back( middlePoint(subdiv[j*3+1], subdiv[j*3+2]) );
                grid.push_back( centroid(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

                triangle2.clear();

                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+2]) );
                triangle2.push_back( middlePoint(subdiv[j*3+1], subdiv[j*3+2]) );
                triangle2.push_back( subdiv[j*3+2] );
                grid.push_back( centroid(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

            }

            subdiv = subdiv2;

        }

    }

    return;

}

GENPoint Radscene::centroid(vector<GENPoint> &shape) {

  double xprime=0.,yprime=0.,zprime=0.;

  for (unsigned int i=0;i<shape.size();i++) {

    xprime += shape[i][0];
    yprime += shape[i][1];
    zprime += shape[i][2];

  }

  xprime/=shape.size();
  yprime/=shape.size();
  zprime/=shape.size();

  return GENPoint::Cartesian(xprime,yprime,zprime);

}

GENPoint Radscene::middlePoint(GENPoint &point1, GENPoint &point2) {

    double xprime, yprime, zprime;

    xprime = 0.5*(point1[0]+point2[0]);
    yprime = 0.5*(point1[1]+point2[1]);
    zprime = 0.5*(point1[2]+point2[2]);

    return GENPoint::Cartesian(xprime, yprime,zprime);

}

void Radscene::averageIrradianceOnSurfaces(vector<unsigned int> count,vector<float> &irrad) {

    for (unsigned int i=0;i<count.size();i++) {

        float average=0.;
        for (unsigned int j=0;j<count[i];j++) {
            average+=irrad.front();
            irrad.erase(irrad.begin());
        }
        average/=float(count[i]);
        irrad.push_back(average);

    }
}

// THE XML SCENE

XmlScene::XmlScene(string inputFile, ostream* pLogFileStr, bool climateFileRequired/*=true*/){

    // logStream is directed by default to the "cout" streambuf
    if(pLogFileStr!=NULL)  // If a logfile stream is provided, redirect logStream to the file stream.
        logStream.rdbuf(pLogFileStr->rdbuf());
    if (!logStream.good())
        throw(string("Unable to define correctly the logStream."));

    // copy the inputFile name
    this->inputFile = inputFile;

    // reads the Radiance file description and put the facades in the surfaceVector
    logStream << "Reading XML file..." << endl << flush;

    TiXmlDocument XMLFile(inputFile.c_str());
    if(!XMLFile.LoadFile()) throw(string("Error loading XML description file. ") + string(" Error #") + toString(XMLFile.ErrorId()) + string(" : ") + toString(XMLFile.ErrorDesc()));

    // browse through the XML file
    TiXmlHandle XMLHandle(&XMLFile);
    TiXmlElement* citysim = XMLHandle.FirstChild("CitySim").ToElement();
    if (!citysim) throw(string("Error in XML file: no CitySim tag."));

    // load the climate
    TiXmlElement* climate = XMLHandle.FirstChild("CitySim").FirstChild("Climate").ToElement();

    if (climateFileRequired && !climate) {
        throw(string("Error in XML file: no Climate tag."));
    }
    else if (climate){
        // Check if climate file exists, first using absolute path, and if it does not exist then with relative path
        string climateFilePath = climate->Attribute("location");
        if (!std::ifstream(climateFilePath.c_str())) {
            // construct path based on .xml file location
            string relativeClimateFilePath = inputFile.substr(0,inputFile.find_last_of("/\\"));
            if (climateFilePath.substr(0,2)=="./")
                climateFilePath = climateFilePath.substr(2,climateFilePath.length()-2); // keep only the name of the file
            while (climateFilePath.substr(0,3)=="../"){
                relativeClimateFilePath = relativeClimateFilePath.substr(0,relativeClimateFilePath.find_last_of("/\\"));
                climateFilePath = climateFilePath.substr(3,climateFilePath.length()-3);
            }
            // gets the delimiter
            string delimiter = (relativeClimateFilePath.find('/') != string::npos) ? "/" : "\\";
            climateFilePath = relativeClimateFilePath + delimiter + climateFilePath;
        }

        // Load climate data only if the file exists
        if(std::ifstream(climateFilePath.c_str())){
            try{
                logStream << "climate file = " << climateFilePath << endl;
                readClimate(climateFilePath);
            }
            catch(string& e){
                logStream << "XmlScene caught exception" << endl;
                if (climateFileRequired) throw e;
                else logStream << "No climate file loaded" << endl;
            }
        }
        else if (climateFileRequired){
            logStream << "ERROR: invalid climate file, the simulation cannot proceed" << endl; // Log as try-catch does not work in release versions
            throw(string("ERROR: invalid climate file, the simulation cannot proceed")); // Throw to stop the program
        }
    }

    // sets the simulation timespan
    unsigned int beginMonth = to<unsigned int>(citysim->FirstChild("Simulation")->ToElement()->Attribute("beginMonth"));
    beginDay   = to<unsigned int>(citysim->FirstChild("Simulation")->ToElement()->Attribute("beginDay"));
    unsigned int endMonth   = to<unsigned int>(citysim->FirstChild("Simulation")->ToElement()->Attribute("endMonth"));
    endDay     = to<unsigned int>(citysim->FirstChild("Simulation")->ToElement()->Attribute("endDay"));
    unsigned int beginYear=0, endYear=0;
    if (citysim->FirstChild("Simulation")->ToElement()->Attribute("beginYear"))
        beginYear  = to<unsigned int>(citysim->FirstChild("Simulation")->ToElement()->Attribute("beginYear"));
    if (citysim->FirstChild("Simulation")->ToElement()->Attribute("endYear"))
        endYear    = to<unsigned int>(citysim->FirstChild("Simulation")->ToElement()->Attribute("endYear"));
    int daysPerMonth[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    // the whole system uses months and days (not years)
    for (unsigned int i = 0; i < beginMonth-1; ++i) beginDay += daysPerMonth[i];
    for (unsigned int i = 0; i < endMonth-1; ++i) endDay += daysPerMonth[i];
    endDay += 365*(endYear-beginYear);

    // sets the seed for the random number generator (must be a uint32_t)
    if ( citysim->FirstChild("Simulation")->ToElement()->Attribute("sprngSeed") ) {
        zigset(to<uint32_t>(citysim->FirstChild("Simulation")->ToElement()->Attribute("sprngSeed")));
        logStream << "sprng seed: " << to<uint32_t>(citysim->FirstChild("Simulation")->ToElement()->Attribute("sprngSeed")) << endl << flush;
    }
    else {
        zigset( static_cast<uint32_t>(26041978) );
        logStream << "sprng seed: " << static_cast<uint32_t>(26041978) << endl << flush;
    }

    // shows the simulation time span
    logStream << "begin day: " << beginDay << "\tend day: " << endDay << endl << flush;

    // creates the district
    pDistrict = new District(XMLHandle.FirstChild("CitySim"), this);

    addAllSurfacesToScene();
}

void XmlScene::addAllSurfacesToScene(){

    // browse the district to create buildings surfaces
    for (unsigned int i=0; i<pDistrict->getnBuildings(); ++i) { // loop on all buildings
        for (unsigned int j=0; j<pDistrict->getBuilding(i)->getnZones(); ++j) { // loop on all zones in the building
            // loop for the walls on this zone
            for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnWalls(); ++k) {
                logStream << "Building " << i << "\tZone: " << j << "\tWall " << k << endl << flush;
                // add the surface to the Buildings surfaces (daylight calculation)
                scene.AddBuildingSurface(GENHandle<Wall>(pDistrict->getBuilding(i)->getZone(j)->getWall(k)));
            }
            // loop for the roofs on this zone
            for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnRoofs(); ++k) {
                logStream << "Building " << i << "\tZone: " << j << "\tRoof " << k << endl << flush;
                // add the surface to the Ground surfaces (meaning NO daylight calculation)
                scene.AddGroundSurface(GENHandle<Roof>(pDistrict->getBuilding(i)->getZone(j)->getRoof(k)));
            }
            // loop on the obstructing surfaces on this zone
            for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnSurfaces(); ++k) {
                logStream << "Building " << i << "\tZone: " << j << "\tSurface " << k << endl << flush;
                // add the surface to the Ground surfaces (meaning NO daylight calculation)
                scene.AddGroundSurface(GENHandle<Surface>(pDistrict->getBuilding(i)->getZone(j)->getSurface(k)));
            }
        }
    }
    logStream << "Buildings' surfaces added to the scene." << endl << flush;

    // browse the district to create ground surfaces
    for (unsigned int i=0; i<pDistrict->getnSurfaces(); ++i) {
        // adds the ground surfaces to the scene
        scene.AddGroundSurface(GENHandle<Surface>(pDistrict->getSurface(i)));
    }
    logStream << "Obstructing surfaces added to the scene." << endl << flush;

    // browse the district to create trees surfaces
    for (size_t i=0; i<pDistrict->getnTrees(); ++i) {
        for (size_t j=0; j<pDistrict->getTree(i)->getnSurfaces(); ++j)
            scene.AddGroundSurface(GENHandle<Surface>(pDistrict->getTree(i)->getSurface(j)));
    }
    logStream << "Trees surfaces added to the scene." << endl;

    // browse the district to create ground surfaces
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        // adds the ground surfaces to the scene
        scene.AddGroundSurface(GENHandle<Ground>(*it));
    }
    logStream << "Ground surfaces added to the scene." << endl;
}

void XmlScene::readClimate(string fileName){

    logStream << "readClimate: " << fileName << endl;
    if(pClimate != NULL){
        delete pClimate;
        pClimate=NULL;
    }

    pClimate = new Climate(fileName,&logStream);
    pClimate->loadCli2(fileName.substr(0,fileName.size()-4) + ".cli2");
    pClimate->loadCli3(fileName.substr(0,fileName.size()-4) + ".cli3");
    climateFile = fileName;

    // shows the information about the location
    logStream << "Location: " << pClimate->getLocation();
    logStream << "\t(Latitude: " << pClimate->getLatitudeN() << " N, Longitude: " << pClimate->getLongitudeE();
    logStream << " E, Altitude: " << pClimate->getAltitude() << " m, Meridian: " << pClimate->getMeridian() << " h)" << endl << flush;

    // gets the information from the climate file
    float latN = pClimate->getLatitudeN();
    float longE = pClimate->getLongitudeE();
    float merE = pClimate->getMeridian()*360./24.;
    float northW = 0.f;

    // initilisation of the location for the sun and the scene
    SKYSiteLocation location(GENAngle::Degrees(latN), GENAngle::Degrees(longE), GENAngle::Degrees(merE), GENAngle::Degrees(northW));
    pSun = new SKYSun(location);
    scene.SetLocation(location);

}

XmlScene::~XmlScene() {

#ifdef DEBUG
    fstream output_IAM("iam_scene.txt", ios::out | ios::trunc);
    output_IAM << ss_IAM.rdbuf();
    output_IAM.close();
    fstream output_IAM_irradiance("iam_scene_irradiance.txt", ios::out | ios::trunc);
    output_IAM_irradiance << ss_IAM_irradiance.rdbuf();
    output_IAM_irradiance.close();
#endif // DEBUG

    //logStream << "Destructor of XmlScene." << endl << flush;
    delete pDistrict;

}

void XmlScene::exportXMLFile(string fileName){
    // define .xml file name
    logStream << "exportXMLFile with fileName=" << fileName << endl;
    if(fileName=="")
        fileName="simulationModel.xml";
    else{
        if(fileName.substr(fileName.size()-4,4)!=".xml"){
            logStream << "exportXML: fileName must end with .xml" << endl << flush;
            throw(string("exportXML: fileName must end with .xml"));
        }
    }
    ofstream file(fileName.c_str());
    if (!file.good()) throw(string("Error creating model .xml file: "+fileName));
    logStream << "file.good()=" << file.good() << endl;

    // define climate file name (same directory as .xml file) and write climate file
    string cliFileName="";
    if(pClimate!=NULL){
        cliFileName=fileName.substr(0,fileName.size()-4)+".cli";
        pClimate->exportCliFile(cliFileName);
    }

    // set the precision for the floating point numbers output in the CitySim XML file
    file.setf(ios::fixed); // set fixed floating format
    file.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    file.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // write .xml file header
    file << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
    file << "<CitySim name=\"test\">" << endl;

    // extract day/month version of beginDay and endDay
    //unsigned int daysPerMonth[] = {31,28,31,30,31,30,31,31,30,31,30,31}; -> const static class attribute
    int bMonth=0,eMonth=0,bDay=beginDay,eDay=endDay;
    while (bDay > daysPerMonth[bMonth]){
        bDay -= daysPerMonth[bMonth];
        ++bMonth;
    }
    while (eDay > daysPerMonth[eMonth]){
        eDay -= daysPerMonth[eMonth];
        ++eMonth;
    }
    file << "\t<Simulation beginMonth=\"" << bMonth+1 << "\" endMonth=\"" << eMonth+1 << "\" beginDay=\"" << bDay << "\" endDay=\"" << eDay << "\"/>" << endl;

    // drop path information to write climate file name in .xml file
    if(cliFileName.find_last_of("/\\")!=std::string::npos)
        cliFileName=cliFileName.substr(cliFileName.find_last_of("/\\")+1);
    file << "\t<Climate location=\""<< cliFileName << "\" city=\"Unknown\"/>" << endl;

    // write scene in .xml file
    pDistrict->writeXML(file, "\t");
    file << "</CitySim>" << endl;
    file.close();

    // write the output results (if they exist) in a binary file
    if (timeStepsSimulated>0) {
        writeSWHeaderText(fileName.substr(0,fileName.size()-4)+"_SW.dat");
        writeSWResultsBinary(fileName.substr(0,fileName.size()-4)+"_SW.dat");
    }

}

void XmlScene::exportGML(string fileName, const vector<double>& origin) {
    // define .xml file name
    //logStream << "exportGMLFile with fileName=" << fileName << endl;
    if(fileName=="") fileName=getInputFileNoExt()+".gml";
    else if(fileName.substr(fileName.size()-4,4)!=".gml") throw(string("exportGML: fileName must end with .gml"));

    ofstream file(fileName.c_str());
    if (!file.good()) throw(string("Error creating model .gml file: "+fileName));

    // set the precision for the floating point numbers output in the CityGML file
    file.setf(ios::fixed); // set fixed floating format
    file.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    file.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // write .xml file header
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    file << "<core:CityModel xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
         << " xsi:schemaLocation=\"http://www.opengis.net/citygml/2.0 http://www.citygmlwiki.org/images/a/ac/EnergyADE.xsd\"\n"
         << " xmlns:xAL=\"urn:oasis:names:tc:ciq:xsdschema:xAL:2.0\"\n"
         << " xmlns:app=\"http://www.opengis.net/citygml/appearance/2.0\" xmlns:wfs=\"http://www.opengis.net/wfs\"\n"
         << " xmlns:genobj=\"http://www.opengis.net/citygml/generics/2.0\"\n"
         << " xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema\"\n"
         << " xmlns:gml=\"http://www.opengis.net/gml\" xmlns:core=\"http://www.opengis.net/citygml/2.0\"\n"
         << " xmlns:group=\"http://www.opengis.net/citygml/cityobjectgroup/2.0\"\n"
         << " xmlns:energy=\"http://www.sig3d.org/citygml/2.0/energy/1.0\"\n"
         << " xmlns:bldg=\"http://www.opengis.net/citygml/building/2.0\">" << endl;
    // write scene in .xml file
    pDistrict->writeGML(file,"\t",origin);
    file << "</core:CityModel>" << endl;
    file.close();
}

void XmlScene::exportDXF(string fileName) {

    if(fileName=="") fileName=inputFile.substr(0,inputFile.size()-4) + ".dxf";
    else if(fileName.substr(fileName.size()-4,4)!=".dxf") fileName.append(".dxf");

    // open the output file
    ofstream outputDxf(fileName.c_str());
    if (!outputDxf.good()) throw(string("Error creating model: " + fileName));

    // setting good precision for the points
    outputDxf.setf(ios::fixed); // set fixed floating format
    outputDxf.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    outputDxf.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // writing the HEADER
    outputDxf << "999\n" << "CitySim" << endl;
    outputDxf << "0\n" << "SECTION" << endl;
    outputDxf << "2\n" << "HEADER" << endl;
    outputDxf << "9\n" << "$ACADVER" << endl;
    outputDxf << "1\n" << "AC1500" << endl; // AutoCAD 2000
    outputDxf << "9\n" << "$INSUNITS" << endl;
    outputDxf << "70\n" << "6" << endl; //6 = Meters
    outputDxf << "0\n" << "ENDSEC" << endl;

    // section of the ENTITIES
    outputDxf << "0\n" << "SECTION\n" << " 2\n" << "ENTITIES" << endl;

    // loop on the grounds -> layer GROUND
    for (forward_list<Ground*>::iterator it=getDistrict()->getGrounds()->begin();it!=getDistrict()->getGrounds()->end();++it) {
        // proceed to the triangulation
        std::shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
        RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
        meshBuilder.AddSurface((*it)->getVertices()->begin(),(*it)->getVertices()->end(),0);
        const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
        for (unsigned int i=0; i<indexedFaceSet->TriangleCount()*3; i+=3) {

            // header of the Triangle in 3DFACE format
            outputDxf << " 0\n" << "3DFACE" << endl;
            outputDxf << "8" << endl; // now comes the layer
            outputDxf << "GROUND" << endl; // layer for the grounds
            outputDxf << "62" << endl; // 62 stands for the color of the surface
            outputDxf << "2" << endl; // yellow for the grounds

            // write those three points and repeat the first one
            for (unsigned int j=0; j<3; ++j) {
                outputDxf << " 1" << j << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+j]))[GENPointCoords::X] << "\n"
                          << " 2" << j << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+j]))[GENPointCoords::Y] << "\n"
                          << " 3" << j << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+j]))[GENPointCoords::Z] << endl;
            }
        }
    } // end of loop on grounds

    // loop on buildings and zones
    for (size_t buildingIndex=0; buildingIndex<pDistrict->getnBuildings(); ++buildingIndex) {
        for (size_t zoneIndex=0; zoneIndex<pDistrict->getBuilding(buildingIndex)->getnZones(); ++zoneIndex) {
            // gets the Surfaces
            vector<Surface*> zoneSurfaces = pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getAllSurfaces();
            for (size_t k=0; k<zoneSurfaces.size();++k) {
                // proceed to the triangulation
                std::shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
                RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
                meshBuilder.AddSurface(zoneSurfaces[k]->getVertices()->begin(),zoneSurfaces[k]->getVertices()->end(),0);
                const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
                for (unsigned int i=0; i<indexedFaceSet->TriangleCount()*3; i+=3) {

                    // header of the Triangle in 3DFACE format
                    outputDxf << " 0\n" << "3DFACE" << endl;
                    outputDxf << "8" << endl; // now comes the layer
                    outputDxf << "BUILDING_" << pDistrict->getBuilding(buildingIndex)->getId() << endl; // one layer for each building starting at 1
                    outputDxf << "62" << endl; // 62 stands for the color of the surface

                    if (zoneSurfaces[k]->getType() == Surface::WALL) outputDxf << "1" << endl; // red for the buildings
                    else if (zoneSurfaces[k]->getType() == Surface::GROUND) outputDxf << "2" << endl; // yellow for the grounds
                    else if (zoneSurfaces[k]->getType() == Surface::ROOF) outputDxf << "3" << endl;
                    else if (zoneSurfaces[k]->getType() == Surface::FLOOR) outputDxf << "4" << endl;
                    else if (zoneSurfaces[k]->getType() == Surface::SURFACE) outputDxf << "5" << endl;

                    // write those three points and repeat the first one
                    for (unsigned int j=0; j<3; ++j) {
                        outputDxf << " 1" << j << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+j]))[GENPointCoords::X] << "\n"
                                  << " 2" << j << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+j]))[GENPointCoords::Y] << "\n"
                                  << " 3" << j << "\n" << GENPoint(indexedFaceSet->GetVertex(pointIndices[i+j]))[GENPointCoords::Z] << endl;
                    }
                }
            } // end of loop on surfaces
        } // end of loop on zones
    } // end of loop on buildings

    outputDxf << " 0\n" << "ENDSEC\n" << " 0\n" << "EOF" << endl;
    outputDxf.close();

}

void XmlScene::exportSTL_ascii(string fileName) {

    if(fileName=="") fileName=inputFile.substr(0,inputFile.size()-4) + ".stl";
    else if(fileName.substr(fileName.size()-4,4)!=".stl") fileName.append(".stl");

    // open the output file
    ofstream output(fileName.c_str());
    if (!output.good()) throw(string("Error creating model: " + fileName));

    // setting good precision for the points
    output.setf(ios::fixed); // set fixed floating format
    output.unsetf(ios::floatfield); //  precision will only specifies the maximum number of digits to be displayed, but not the minimum
    output.precision(numeric_limits<float>::max_digits10); // set the precision to the maximum digits possible with float

    // loop on buildings and zones
    for (size_t buildingIndex=0; buildingIndex<pDistrict->getnBuildings(); ++buildingIndex) {
        for (size_t zoneIndex=0; zoneIndex<pDistrict->getBuilding(buildingIndex)->getnZones(); ++zoneIndex) {
            // gets the Surfaces
            vector<Surface*> zoneSurfaces = pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getAllSurfaces();
            for (size_t k=0; k<zoneSurfaces.size();++k) {
                // proceed to the triangulation
                shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
                RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
                meshBuilder.AddSurface(zoneSurfaces[k]->getVertices()->begin(),zoneSurfaces[k]->getVertices()->end(),0);
                const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
                // outputs the name of the surface
                output << "solid b" << pDistrict->getBuilding(buildingIndex)->getId() << "_s" << zoneSurfaces[k]->getId() << endl;
                for (unsigned int i=0; i<indexedFaceSet->TriangleCount(); ++i) {
                    // computes the normal
                    GENPoint triangleNormal = GENPoint::UnitVector(GEOMPolygonInfo::NormalFromThreePoints(GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3  ])),
                                                                                                          GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+1])),
                                                                                                          GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+2]))));

                    // each triangle is named a facet and has a normal
                    output << "facet normal " << triangleNormal[GENPointCoords::X] << " "
                                              << triangleNormal[GENPointCoords::Y] << " "
                                              << triangleNormal[GENPointCoords::Z] << " " << endl;

                    output << "\touter loop" << endl;
                    // write those three points and repeat the first one
                    for (unsigned int j=0; j<3; ++j) {
                        output << "\t\tvertex " << GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+j]))[GENPointCoords::X] << " "
                                                << GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+j]))[GENPointCoords::Y] << " "
                                                << GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+j]))[GENPointCoords::Z] << endl;
                    }
                    output << "\tendloop" << endl;
                    output << "endfacet" << endl;
                }
                output << "endsolid b" << pDistrict->getBuilding(buildingIndex)->getId() << "_s" << zoneSurfaces[k]->getId() << endl;
            } // end of loop on surfaces
        } // end of loop on zones
    } // end of loop on buildings

    // loop on the grounds
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        // proceed to the triangulation
        shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
        RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
        meshBuilder.AddSurface((*it)->getVertices()->begin(),(*it)->getVertices()->end(),0);
        const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
        // outputs the name of the surface
        output << "solid g" << (*it)->getId() << endl;
        for (unsigned int i=0; i<indexedFaceSet->TriangleCount(); ++i) {
            // computes the normal
            GENPoint triangleNormal = GENPoint::UnitVector(GEOMPolygonInfo::NormalFromThreePoints(GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3  ])),
                                                                                                  GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+1])),
                                                                                                  GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+2]))));

            // each triangle is named a facet and has a normal
            output << "facet normal " << triangleNormal[GENPointCoords::X] << " "
                                      << triangleNormal[GENPointCoords::Y] << " "
                                      << triangleNormal[GENPointCoords::Z] << " " << endl;

            output << "\touter loop" << endl;
            // write those three points and repeat the first one
            for (unsigned int j=0; j<3; ++j) {
                output << "\t\tvertex " << GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+j]))[GENPointCoords::X] << " "
                                        << GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+j]))[GENPointCoords::Y] << " "
                                        << GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+j]))[GENPointCoords::Z] << endl;
            }
            output << "\tendloop" << endl;
            output << "endfacet" << endl;
        }
        output << "endsolid g" << (*it)->getId() << endl;
    } // end of loop on surfaces

    output.close();

}

void XmlScene::exportSTL_binary(string fileName) {

//    UINT8[80] - Header
//    UINT32 - Number of triangles

//    foreach triangle
//    REAL32[3] - Normal vector
//    REAL32[3] - Vertex 1
//    REAL32[3] - Vertex 2
//    REAL32[3] - Vertex 3
//    UINT16 - Attribute byte count
//    end

    if(fileName=="") fileName=inputFile.substr(0,inputFile.size()-4) + ".stl";
    else if(fileName.substr(fileName.size()-4,4)!=".stl") fileName.append(".stl");

    // open the stringstream for writing
    ostringstream output;

    // loop on all triangles
    shared_ptr<RENIndexedFaceSet> pIndexedFaceSet = scene.IndexedFaceSet();
    // get bounding Z
    pair<float,float> Zmin_max = pIndexedFaceSet->getBoundingZ();
    // get Indexed Triangulated Surfaces
    const std::vector<unsigned int>& pointIndices=pIndexedFaceSet->GetIndices();
    // initialize the surfaceCount and a float value to be used for output
    uint32_t surfaceCount = 0;
    float value;
    for (unsigned int i=0; i<pIndexedFaceSet->TriangleCount(); ++i) {

        // computes the normal
        GENPoint triangleNormal = GENPoint::UnitVector(GEOMPolygonInfo::NormalFromThreePoints(GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3  ])),
                                                                                              GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3+1])),
                                                                                              GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3+2]))));

        // write the normal of each triangle
        value = triangleNormal[GENPointCoords::X];
        output.write(reinterpret_cast<char*>(&value),sizeof(value));
        value = triangleNormal[GENPointCoords::Y];
        output.write(reinterpret_cast<char*>(&value),sizeof(value));
        value = triangleNormal[GENPointCoords::Z];
        output.write(reinterpret_cast<char*>(&value),sizeof(value));

        // write the vertices
        for (unsigned int k=0; k<3; k++) {
            value = GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3]))[k];
            output.write(reinterpret_cast<char*>(&value),sizeof(value));
        }
        for (unsigned int k=0; k<3; k++) {
            value = GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3+1]))[k];
            output.write(reinterpret_cast<char*>(&value),sizeof(value));
        }
        for (unsigned int k=0; k<3; k++) {
            value = GENPoint(pIndexedFaceSet->GetVertex(pointIndices[i*3+2]))[k];
            output.write(reinterpret_cast<char*>(&value),sizeof(value));
        }

        // write the Attribute byte count (the id of the triangle)
        output.write(reinterpret_cast<char*>(&i),sizeof(uint16_t));

        // increment the surfaceCount
        ++surfaceCount;
    }

    // loop on the floors at 0 m
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the floors
            for (vector<Floor*>::iterator it=pDistrict->getBuilding(j)->getZone(zone)->getFloors()->begin();it!=pDistrict->getBuilding(j)->getZone(zone)->getFloors()->end();++it) {
                // proceed to the triangulation
                shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
                RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
                meshBuilder.AddSurface((*it)->getVertices()->begin(),(*it)->getVertices()->end(),0);
                const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
                // outputs the triangles of the ground
                for (unsigned int i=0; i<indexedFaceSet->TriangleCount(); ++i) {

                    // computes the normal, reverse the vertices
                    GENPoint triangleNormal = GENPoint::UnitVector(GEOMPolygonInfo::NormalFromThreePoints(GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+2])),
                                                                                                          GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+1])),
                                                                                                          GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3  ]))));

                    // write the normal of each triangle
                    value = triangleNormal[GENPointCoords::X];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                    value = triangleNormal[GENPointCoords::Y];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                    value = triangleNormal[GENPointCoords::Z];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));

                    // write the vertices
                    for (unsigned int k=0; k<3; k++) {
                        value = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+2]))[k];
                        // deduce 1/5 of the maximal height of the buildings to the min value
                        if (k==2) value = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);
                        output.write(reinterpret_cast<char*>(&value),sizeof(value));
                    }
                    for (unsigned int k=0; k<3; k++) {
                        value = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+1]))[k];
                        // deduce 1/5 of the maximal height of the buildings to the min value
                        if (k==2) value = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);
                        output.write(reinterpret_cast<char*>(&value),sizeof(value));
                    }
                    for (unsigned int k=0; k<3; k++) {
                        value = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3  ]))[k];
                        // deduce 1/5 of the maximal height of the buildings to the min value
                        if (k==2) value = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);
                        output.write(reinterpret_cast<char*>(&value),sizeof(value));
                    }
                    // write the Attribute byte count (the id of the triangle)
                    value = 0.f;
                    output.write(reinterpret_cast<char*>(&value),sizeof(uint16_t));

                    // increment the surfaceCount
                    ++surfaceCount;
                }
            }
        }
    }

    // loop on the grounds at 0 m
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        // proceed to the triangulation
        shared_ptr<RENIndexedFaceSet> indexedFaceSet(new RENIndexedFaceSet());
        RENIndexedFaceSetBuilder meshBuilder(indexedFaceSet);
        meshBuilder.AddSurface((*it)->getVertices()->begin(),(*it)->getVertices()->end(),0);
        const std::vector<unsigned int>& pointIndices=indexedFaceSet->GetIndices();
        // outputs the triangles of the ground
        for (unsigned int i=0; i<indexedFaceSet->TriangleCount(); ++i) {

            // computes the normal, reverse the vertices
            GENPoint triangleNormal = GENPoint::UnitVector(GEOMPolygonInfo::NormalFromThreePoints(GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+2])),
                                                                                                  GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+1])),
                                                                                                  GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3  ]))));

            // write the normal of each triangle
            value = triangleNormal[GENPointCoords::X];
            output.write(reinterpret_cast<char*>(&value),sizeof(value));
            value = triangleNormal[GENPointCoords::Y];
            output.write(reinterpret_cast<char*>(&value),sizeof(value));
            value = triangleNormal[GENPointCoords::Z];
            output.write(reinterpret_cast<char*>(&value),sizeof(value));

            // write the vertices
            for (unsigned int k=0; k<3; k++) {
                value = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+2]))[k];
                // deduce 1/5 of the maximal height of the buildings to the min value
                if (k==2) value = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);
                output.write(reinterpret_cast<char*>(&value),sizeof(value));
            }
            for (unsigned int k=0; k<3; k++) {
                value = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+1]))[k];
                // deduce 1/5 of the maximal height of the buildings to the min value
                if (k==2) value = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);
                output.write(reinterpret_cast<char*>(&value),sizeof(value));
            }
            for (unsigned int k=0; k<3; k++) {
                value = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3  ]))[k];
                // deduce 1/5 of the maximal height of the buildings to the min value
                if (k==2) value = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);
                output.write(reinterpret_cast<char*>(&value),sizeof(value));
            }
            // write the Attribute byte count (the id of the triangle)
            value = 0.f;
            output.write(reinterpret_cast<char*>(&value),sizeof(uint16_t));

            // increment the surfaceCount
            ++surfaceCount;
        }
        // outputs the prisms of the ground
        for (unsigned int i=0; i<indexedFaceSet->TriangleCount(); ++i) {

            // loop on the 3 segments of the triangles
            for (unsigned int j=0; j<3; ++j) {

                // A' and B' are respectively the 3 segments of the triangles
                GENPoint Aprime = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+(j)]));
                GENPoint Bprime = GENPoint(indexedFaceSet->GetVertex(pointIndices[i*3+(j+1)%3]));

                // first triangle; A', B'(0 m), B'
                GENPoint Bprime0 = Bprime;
                // deduce 1/5 of the maximal height of the buildings to the min value
                Bprime0[GENPointCoords::Z] = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);

                // computes the normal
                GENPoint triangleNormal = GENPoint::UnitVector(GEOMPolygonInfo::NormalFromThreePoints(Aprime, Bprime0, Bprime));

                // write the normal of each triangle
                value = triangleNormal[GENPointCoords::X];
                output.write(reinterpret_cast<char*>(&value),sizeof(value));
                value = triangleNormal[GENPointCoords::Y];
                output.write(reinterpret_cast<char*>(&value),sizeof(value));
                value = triangleNormal[GENPointCoords::Z];
                output.write(reinterpret_cast<char*>(&value),sizeof(value));

                // write the vertices
                for (unsigned int k=0; k<3; k++) {
                    value = Aprime[k];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                }
                for (unsigned int k=0; k<3; k++) {
                    value = Bprime0[k];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                }
                for (unsigned int k=0; k<3; k++) {
                    value = Bprime[k];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                }

                // write the Attribute byte count (the id of the triangle)
                value = 0.f;
                output.write(reinterpret_cast<char*>(&value),sizeof(uint16_t));

                // second triangle; A', A'(0 m), B'(0 m)
                GENPoint Aprime0 = Aprime;
                // deduce 1/5 of the maximal height of the buildings to the min value
                Aprime0[GENPointCoords::Z] = Zmin_max.first - 0.2*(Zmin_max.second-Zmin_max.first);

                // computes the normal
                triangleNormal = GENPoint::UnitVector(GEOMPolygonInfo::NormalFromThreePoints(Aprime, Aprime0, Bprime0));

                // write the normal of each triangle
                value = triangleNormal[GENPointCoords::X];
                output.write(reinterpret_cast<char*>(&value),sizeof(value));
                value = triangleNormal[GENPointCoords::Y];
                output.write(reinterpret_cast<char*>(&value),sizeof(value));
                value = triangleNormal[GENPointCoords::Z];
                output.write(reinterpret_cast<char*>(&value),sizeof(value));

                // write the vertices
                for (unsigned int k=0; k<3; k++) {
                    value = Aprime[k];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                }
                for (unsigned int k=0; k<3; k++) {
                    value = Aprime0[k];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                }
                for (unsigned int k=0; k<3; k++) {
                    value = Bprime0[k];
                    output.write(reinterpret_cast<char*>(&value),sizeof(value));
                }

                // write the Attribute byte count (the id of the triangle)
                value = 0.f;
                output.write(reinterpret_cast<char*>(&value),sizeof(uint16_t));

                // increment the surfaceCount
                surfaceCount+=2; // two triangles
            }
        }
    } // end of loop on surfaces

    // open the output file
    ofstream outputFile(fileName.c_str(),ios::binary | ios::trunc);
    if (!outputFile.good()) throw(string("Error creating model: " + fileName));

    // write the header
    char header[80] = "Binary STL from CitySim";
    outputFile.write(header,sizeof(header));

    // write the number of triangles
    outputFile.write(reinterpret_cast<char*>(&surfaceCount),sizeof(surfaceCount));

    // write the stringstream and close the file
    outputFile << output.str();
    outputFile.close();

}

void XmlScene::exportRAD(string fileName) {

    exportRadFile(fileName);
    exportInpFile(fileName);

}

void XmlScene::initialiseFarField() {

    // check the size of the far field obstructions vector, no obstructions -> return without changes in lv
    if (pDistrict->getnFarFieldObstructions() == 0) return;

    // creates the vector containing the Far Field obstructions
    farFieldOccludedPatchFraction.assign(tregenzaSky.PatchCount()/2, 0.);

    // create the used variables
    double phi1,phi2,theta1,theta2;
    double x1,y1,x2,y2; // this comes from the XML file (as float)
    double phiInter1,phiInter2,thetaInter1,thetaInter2;
    double phiStart,thetaStart,phiStop,thetaStop;

    // vector containing the occluded patch fraction
    //ostringstream outputFarField;
    //outputFarField << "#Patch number\tOccluded Spherical Area\tZone Spherical Area\tPatch Spherical Area\tOccludedPatchFraction" << endl;

    // for debugging reasons
    //ostringstream outputScreen;

    // loop on the bands
    for (unsigned int band=0;band<tregenzaSky.getBands();++band) {
        // loop on the zones per band
        for (unsigned int zone=0;zone<tregenzaSky.getPatchesPerBand(band)+1;++zone) {

            // get the corresponding patch number
            int patchNumber = tregenzaSky.getPatchIndex(band,zone%tregenzaSky.getPatchesPerBand(band));

            // computes the phi and theta according to the patch number
            phi1=max((360.*(static_cast<double>(zone)-0.5))/tregenzaSky.getPatchesPerBand(band),0.);
            phi2=min((360.*(static_cast<double>(zone)+0.5))/tregenzaSky.getPatchesPerBand(band),360.);
            theta1=band*12.;
            theta2=min((band+1)*12., 90.);

            // if phi1=phi2, then not a zone
            if (phi1==phi2) break;

            // some verifications
            //outputScreen << "band: " << band << " zone: " << zone << " patch number: " << patchNumber << endl;
            //outputScreen << " phi1: " << phi1 << " phi2: " << phi2 << " theta1: " << theta1 << " theta2: " << theta2 << endl;

            // loop on the segments for this patch to get the occluded area
            double occludedPatchSphericalArea=0.;

            for (unsigned int segment=0;segment+1<pDistrict->getnFarFieldObstructions();segment++) {

                // get the start and stop point of the corresponding segment (x1,y1) and (x2,y2)
                x1=pDistrict->getFarFieldObstructions(segment).first;
                y1=pDistrict->getFarFieldObstructions(segment).second;
                x2=pDistrict->getFarFieldObstructions(segment+1).first;
                y2=pDistrict->getFarFieldObstructions(segment+1).second;

                // test if this segment is relevant for the current patch, if not then continue to the next segment
                if ( !(x1<=phi2 && x2 >= phi1) ) continue;

                // compute the intersection with phi = constant
                thetaInter1 = (y2-y1)/(x2-x1)*phi1+(y1*x2-y2*x1)/(x2-x1);
                thetaInter2 = (y2-y1)/(x2-x1)*phi2+(y1*x2-y2*x1)/(x2-x1);

                // case 1) when the slope of the segment is positive
                if ( y2>y1 ) {

                    //outputScreen << "Positive slope" << endl;

                    // compute the intersection with theta1, theta2 lines
                    phiInter1 = (theta1-(y1*x2-y2*x1)/(x2-x1))*(x2 -x1)/(y2-y1);
                    phiInter2 = (theta2-(y1*x2-y2*x1)/(x2-x1))*(x2 -x1)/(y2-y1);

                    //outputScreen << "inter1(" << phiInter1 << "," << thetaInter1 << ") - inter2(" << phiInter2 << "," << thetaInter2 << ")" << endl;

                    // verifies if this segment is above the current patch, then use the whole patch area
                    if ( x1 <= phi2 && x2 >= phi1 && thetaInter1 >= theta2 && thetaInter2 >= theta2 ) {
                        occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(max(x1,phi1)/180.*M_PI,theta1/180.*M_PI,min(x2,phi2)/180.*M_PI,theta2/180.*M_PI);
                        //outputScreen << "Whole patch added." << endl;
                        continue;
                    }

                    // compute the intersection with the current patch
                    phiStart   = min(max(min(max(phi1,phiInter1),phi2),x1),x2);
                    thetaStart = min(max(min(max(theta1,thetaInter1),theta2),y1),y2);
                    phiStop    = min(max(max(min(phi2,phiInter2),phi1),x1),x2);
                    thetaStop  = min(max(max(min(theta2,thetaInter2),theta1),y1),y2);

                    //outputScreen << "start(" << phiStart << "," << thetaStart << ") - stop(" << phiStop << "," << thetaStop << ")" << endl;

                    // remove the case where start=stop
                    if (phiStart == phiStop && thetaStart == thetaStop) continue;

                    // compute the occluded patch area
                    occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(phiStart/180.*M_PI,theta1/180.*M_PI,phiStop/180.*M_PI,thetaStart/180.*M_PI)
                                                  +tregenzaSky.triangleSphericalArea(phiStart/180.*M_PI,thetaStart/180.*M_PI,phiStop/180.*M_PI,thetaStop/180.*M_PI);

                    // if the line cuts theta2, then add the rest of the patch
                    if (phiStop != x2) {
                        occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(phiStop/180.*M_PI,theta1/180.*M_PI,phi2/180.*M_PI,theta2/180.*M_PI);
                    }

                    //outputScreen << "Occluded patch area: " << occludedPatchSphericalArea << endl;

                }
                else if ( y2<y1 ) {

                    //outputScreen << "Negative slope" << endl;

                    // compute the intersection with theta1, theta2 lines
                    phiInter1 = (theta2-(y1*x2-y2*x1)/(x2-x1))*(x2 -x1)/(y2-y1); // inversion of theta1 and theta2
                    phiInter2 = (theta1-(y1*x2-y2*x1)/(x2-x1))*(x2 -x1)/(y2-y1); // inversion of theta1 and theta2

                    //outputScreen << "inter1(" << phiInter1 << "," << thetaInter1 << ") - inter2(" << phiInter2 << "," << thetaInter2 << ")" << endl;

                    // verifies if this segment is above the current patch, then use the whole patch area
                    if ( x1 <= phi2 && x2 >= phi1 && thetaInter1 >= theta2 && thetaInter2 >= theta2 ) {
                        occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(max(x1,phi1)/180.*M_PI,theta1/180.*M_PI,min(x2,phi2)/180.*M_PI,theta2/180.*M_PI);
                        continue;
                    }

                    // compute the intersection with the current patch
                    phiStart   = min(max(min(max(phi1,phiInter1),phi2),x1),x2);
                    thetaStart = max(min(min(max(theta1,thetaInter1),theta2),y1),y2); // inversion of the first two functions
                    phiStop    = min(max(max(min(phi2,phiInter2),phi1),x1),x2);
                    thetaStop  = max(min(max(min(theta2,thetaInter2),theta1),y1),y2); // inversion of the first two functions

                    //outputScreen << "start(" << phiStart << "," << thetaStart << ") - stop(" << phiStop << "," << thetaStop << ")" << endl;

                    // remove the case where start=stop
                    if (phiStart == phiStop && thetaStart == thetaStop) continue;

                    // compute the occluded patch area
                    occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(phiStart/180.*M_PI,theta1/180.*M_PI,phiStop/180.*M_PI,thetaStop/180.*M_PI)
                                                  +tregenzaSky.triangleSphericalArea(phiStart/180.*M_PI,thetaStart/180.*M_PI,phiStop/180.*M_PI,thetaStop/180.*M_PI);

                    // if the cut is in theta2, then add a part
                    if (phiStart != x1) {
                        occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(phi1/180.*M_PI,theta1/180.*M_PI,phiStart/180.*M_PI,theta2/180.*M_PI);
                    }

                    //outputScreen << "Occluded patch area: " << occludedPatchSphericalArea << endl;

                }
                else { // y1 and y2 are equal, flat line or segment

                    //outputScreen << "Flat line" << endl;

                    // verifies if this segment is above the current patch, then use the whole patch area
                    if ( x1 <= phi2 && x2 >= phi1 && thetaInter1 >= theta2 && thetaInter2 >= theta2 ) {
                        occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(max(x1,phi1)/180.*M_PI,theta1/180.*M_PI,min(x2,phi2)/180.*M_PI,theta2/180.*M_PI);
                        continue;
                    }

                    // compute the intersection with the patch
                    phiStart   = min(max(phi1,x1),x2);
                    phiStop    = min(max(phi2,x1),x2);
                    thetaStop  = max(min(theta2,thetaInter2),theta1);

                    //outputScreen << "start(" << phiStart << "," << thetaStop << ") - stop(" << phiStop << "," << thetaStop << ")" << endl;

                    // compute the occluded patch area
                    if (thetaStop != theta1) {
                        occludedPatchSphericalArea += tregenzaSky.zoneSphericalArea(phiStart/180.*M_PI,theta1/180.*M_PI,phiStop/180.*M_PI,thetaStop/180.*M_PI);
                    }
                    //outputScreen << "Occluded patch area: " << occludedPatchSphericalArea << endl;

                }

            } // end the loop on the segments

            // compute the fraction of the spherical area, in two parts for the extreme patches
            farFieldOccludedPatchFraction[patchNumber] += occludedPatchSphericalArea / tregenzaSky.GetPatch(patchNumber)->solidAngle();

            // saves the computation for verification
            //outputFarField << patchNumber << "\t" << occludedPatchSphericalArea << "\t" << tregenzaSky.zoneSphericalArea(phi1/180.*M_PI,theta1/180.*M_PI,phi2/180.*M_PI,theta2/180.*M_PI)
            //               << "\t" << tregenzaSky.GetPatch(patchNumber)->solidAngle() << "\t" << farFieldOccludedPatchFraction[patchNumber] << endl;

        } // end of the loop on the patches of the band
    } // end of the loop on bands

    // debugging outputs
    //save(string("farField.txt"),outputFarField);
    //save(string("farFieldDetails.txt"),outputScreen);
    //save("farFieldOcclPatchFraction.txt",farFieldOccludedPatchFraction);
    //save("farFieldObstructions.txt",pDistrict->getFarFieldObstructions());

    return;
}

void XmlScene::computeFarField() {

    // check the size of the far field obstructions vector, no obstructions -> return without changes in lv
    if (pDistrict->getnFarFieldObstructions() == 0) return;

    // loop on the Tregenza zones in the sky
    for (unsigned int patchNumber=0; patchNumber<tregenzaSky.PatchCount()/2; ++patchNumber) {
        // weight the cumulative radiance between the sky and ground
        lv[patchNumber] = lv[patchNumber]*(1.-farFieldOccludedPatchFraction[patchNumber]) + groundRadiance*farFieldOccludedPatchFraction[patchNumber];
    }

    return;
}

bool XmlScene::sunVisibleFarField(float sunAzimuth, float sunElevation) {

    // checks if the azimuth is in the correct range to be compatible with Far Field data
    sunAzimuth = fmod(sunAzimuth+360.f,360.f);

    float x1,y1,x2,y2;

    // checks if the sun is visible according to the Far Field obstructions
    for (unsigned int segment=0;segment+1<pDistrict->getnFarFieldObstructions();segment++) {

        // get the start and stop point of the corresponding segment (x1,y1) and (x2,y2)
        x1=pDistrict->getFarFieldObstructions(segment).first;
        x2=pDistrict->getFarFieldObstructions(segment+1).first;

        // for each segment checks if the sunAzimuth and sunElevation is higher
        if ( sunAzimuth >= x1 && sunAzimuth < x2 ) { // we are in

            y1=pDistrict->getFarFieldObstructions(segment).second;
            y2=pDistrict->getFarFieldObstructions(segment+1).second;

            // computes the elevation corresponding to the sun azimuth
            float y3 = (y2-y1)/(x2-x1)*sunAzimuth+(y1*x2-y2*x1)/(x2-x1);

            return (sunElevation > y3);
        }
    }
    return true;
}

void XmlScene::computeShortWave_Beam(unsigned int day, unsigned int hour) {

    // initialise the beam vector
    irradiationBeam.assign(scene.SurfaceCount(), 0.f);

    // initialisation of the sun (for the VFC)
    pSun->SetDay(day);
    if (pSun->SetClockTime1(static_cast<float>(hour))) { // if the day has started

        // gets Ibn from the climate file
        float Ibn = pClimate->getIbn(day,hour);

        // computes the Igh
        float Igh = 0.f;
        // direct part contribution to Igh, if sun is above horizon
        if (sunVisibleFarField(pSun->GetPosition().Azimuth().degrees(), pSun->GetPosition().Altitude().degrees())) {
            float cosTheta = GEN::dot_product(GENPoint::Cartesian(0.f,0.f,1.f),pSun->GetPosition());
            if (cosTheta > 0.f) {
                Igh += Ibn * cosTheta;
            }
        }
        pClimate->addIgh(Igh);

        // now the main loop on surfaces
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
            // direct part
            if (sunVisibleFarField(pSun->GetPosition().Azimuth().degrees(), pSun->GetPosition().Altitude().degrees()) && Ibn > 0.f) {
                double cosTheta = GEN::dot_product(scene.GetSurface(surfaceIndex).Normal(),pSun->GetPosition());
                if (cosTheta > 0.) {

                    irradiationBeam[surfaceIndex] = scene.GetSurface(surfaceIndex).InsolationFactors().GetInsolationFactor(*pSun) * Ibn * cosTheta; // fraction of the surface that is light by the sun

                }
            }
        } // end the loop on the surfaces
    } // end if sun is up

    // end of the calculations, save in all surfaces the irradiation
    for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
    {
        // saves all the irradiance
        ((Surface*)(it->SurfaceDelegate()))->addShortWaveIrradiance(irradiationBeam[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

        // saves all the beam irradiance
        ((Surface*)(it->SurfaceDelegate()))->setBeamIrradiance(irradiationBeam[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

        // saves the angle between the sun and the surface
        ((Surface*)(it->SurfaceDelegate()))->setBeamAngle(acosf(GEN::dot_product(it->Normal(),pSun->GetPosition())/(it->Normal().Radius()*pSun->GetPosition().Radius()))/M_PI*180.f);

    }

    return;

}

void XmlScene::computeShortWave_Diffuse() {

    // initialise the diffuse vectors
    irradiationDiffuseSky.assign(scene.SurfaceCount(), 0.f);
    irradiationDiffuseGround.assign(scene.SurfaceCount(), 0.f);

    // we assume the computeRadiance and computeFarField are launched beforehand
    // lv[i] and groundRadiance

    // computes the contribution to Igh
    float Igh = 0.f;
    // diffuse part
    for (unsigned int patchNo = 0; patchNo < tregenzaSky.PatchCount()/2; ++patchNo) // loop on the sky patches
    {
        Igh += tregenzaSky.GetPatch(patchNo)->formFactor(GENPoint::Cartesian(0.f,0.f,1.f)) * lv[patchNo];
    }
    pClimate->addIgh(Igh);

    // now the main loop on surfaces
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
        // diffuse part sky
        for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).SWViewFactors().GetVFs();
            factors!=scene.GetSurface(surfaceIndex).SWViewFactors().GetLastVF();
            ++factors) // loop on the patches that are non-zero
        {
            // the factors->unobstructed contains the patchSolidAngle and the cos of the angles between the surface normal and the center of patch
            if (factors->unobstructed > 0.f) {
                if ( factors->patchNo < static_cast<int>(tregenzaSky.PatchCount()/2) ) {
                    // diffuse part sky
                    irradiationDiffuseSky[surfaceIndex] += factors->unobstructed * lv[factors->patchNo];
                }
            }
        }
        // diffuse part ground
        irradiationDiffuseGround[surfaceIndex] += ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_ground() * groundRadiance;
    } // end the loop on the surfaces

    // end of the calculations, save in all surfaces the irradiation
    for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
    {
        // saves all the irradiance
        ((Surface*)(it->SurfaceDelegate()))->addShortWaveIrradiance(irradiationDiffuseSky[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]
                                                                    +irradiationDiffuseGround[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

    }

    return;

}

void XmlScene::computeShortWave_Interreflected() {

    // irradiationSW: irradiation without reflections, irradiationSWn: irradiation with reflections
    vector<float> irradiationSW(scene.SurfaceCount(), 0.f), irradiationSWn;

    // end of the calculations, save in all surfaces the irradiation
    for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
    {
        // saves all the irradiance
        irradiationSW[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]=((Surface*)(it->SurfaceDelegate()))->getShortWaveIrradiance();

    }

    // adds the reflections on the surfaces
    irradiationSWn = irradiationSW;
    for (unsigned int r = 0; r < mNbReflections; ++r) {

        // using the CRS sparse matrix format, multiplication of the matrix by the former vector
        vector<float> irradiationSWnew(irradiationSW);
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int i=0; i < getnAi(); ++i) { // loop on the number of elements
            for (unsigned int index=getAi(i); index < getAi(i+1); ++index) {
                irradiationSWnew[i] += getAn(index) * irradiationSWn[getAj(index)];
            }
        }
        // saves the new irradiation vector in the current
        irradiationSWn = irradiationSWnew;

    }

    // end of the calculations, save in all surfaces the irradiation
    for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
    {
        // saves all the irradiance
        ((Surface*)(it->SurfaceDelegate()))->addShortWaveIrradiance(irradiationSWn[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]
                                                                    -irradiationSW[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

    }

    return;

}

void XmlScene::computeShortWave(unsigned int day, unsigned int hour) {

    // irradiationSW: irradiation without reflections, irradiationSWn: irradiation with reflections
    vector<float> irradiationSW(scene.SurfaceCount(), 0.f), irradiationSWn(scene.SurfaceCount(), 0.f);
    vector<float> irradiationSW_IAM(scene.SurfaceCount(), 0.f);
    irradiationBeam.assign(scene.SurfaceCount(), 0.f);
    irradiationDiffuseSky.assign(scene.SurfaceCount(), 0.f);
    irradiationDiffuseGround.assign(scene.SurfaceCount(), 0.f);

    // initialisation of the sun (for the VFC)
    pSun->SetDay(day);
    if (pSun->SetClockTime1(static_cast<float>(hour))) { // if the day has started

        // gets Idh and Ibn from the climate file
        float Idh = pClimate->getIdh(day,hour);
        float Ibn = pClimate->getIbn(day,hour);

        // computes the patches radiance without and then with obstructions
        computeRadiance(day,Idh,Ibn,pDistrict->getGroundAlbedo());
        //save("cumRadiance_" + toString(day) + "_" + toString(hour) + ".txt", lv);
        computeFarField();
        //save("cumRadianceFarField_" + toString(day) + "_" + toString(hour) + ".txt", lv);

        // computes the Igh
        float Igh = 0.f;
        // direct part contribution to Igh, if sun is above horizon
        if (sunVisibleFarField(pSun->GetPosition().Azimuth().degrees(), pSun->GetPosition().Altitude().degrees())) {
            float cosTheta = GEN::dot_product(GENPoint::Cartesian(0.f,0.f,1.f),pSun->GetPosition());
            if (cosTheta > 0.f) {
                Igh += pClimate->getIbn(day,hour) * cosTheta;
            }
        }
        // diffuse part
        for (unsigned int patchNo = 0; patchNo < tregenzaSky.PatchCount()/2; ++patchNo) // loop on the sky patches
        {
            Igh += tregenzaSky.GetPatch(patchNo)->formFactor(GENPoint::Cartesian(0.f,0.f,1.f)) * lv[patchNo];
        }
        pClimate->setIgh(Igh);
        // end of Igh compute

        // now the main loop on surfaces
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
            // diffuse part sky
            for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).SWViewFactors().GetVFs();
                factors!=scene.GetSurface(surfaceIndex).SWViewFactors().GetLastVF();
                ++factors) // loop on the patches that are non-zero
            {
                // the factors->unobstructed contains the patchSolidAngle and the cos of the angles between the surface normal and the center of patch
                if (factors->unobstructed > 0.f) {
                    if ( factors->patchNo < static_cast<int>(tregenzaSky.PatchCount()/2) ) {
                        // diffuse part sky
                        irradiationSW[surfaceIndex] += factors->unobstructed * lv[factors->patchNo];
                        irradiationDiffuseSky[surfaceIndex] += factors->unobstructed * lv[factors->patchNo];
#ifdef DEBUG
                        if (((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()) {
                            float theta = acos(cosAngleBetween(scene.GetSurface(surfaceIndex).Normal(),tregenzaSky.GetPatch(factors->patchNo)->averageVisibleDirection(scene.GetSurface(surfaceIndex).Normal())));
                            ss_IAM << day << "\t" << hour << "\t" << surfaceIndex << "\t"
                                   << scene.GetSurface(surfaceIndex).Normal().Altitude().degrees() << "\t"
                                   << scene.GetSurface(surfaceIndex).Normal().Azimuth().degrees() << "\t"
                                   << "SKY\t" << factors->patchNo + 1 << "\t"
                                   << tregenzaSky.GetPatch(factors->patchNo)->averageVisibleDirection(scene.GetSurface(surfaceIndex).Normal()).Altitude().degrees() << "\t"
                                   << tregenzaSky.GetPatch(factors->patchNo)->averageVisibleDirection(scene.GetSurface(surfaceIndex).Normal()).Azimuth().degrees() << "\t"
                                   << theta*180./M_PI << "\t"
                                   << factors->unobstructed * lv[factors->patchNo] << "\t"
                                   << ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()->getIAM(theta) << endl;
                            irradiationSW_IAM[surfaceIndex] += factors->unobstructed * lv[factors->patchNo] * ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()->getIAM(theta);
                        }
#endif
                    }
#ifdef DEBUG
                    else if (((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()) {
                        // diffuse part ground
                        float theta = acos(cosAngleBetween(scene.GetSurface(surfaceIndex).Normal(),tregenzaSky.GetPatch(factors->patchNo)->averageVisibleDirection(scene.GetSurface(surfaceIndex).Normal())));
                        ss_IAM << day << "\t" << hour << "\t" << surfaceIndex << "\t"
                        << scene.GetSurface(surfaceIndex).Normal().Altitude().degrees() << "\t"
                        << scene.GetSurface(surfaceIndex).Normal().Azimuth().degrees() << "\t"
                        << "GROUND\t" << factors->patchNo + 1 << "\t"
                        << tregenzaSky.GetPatch(factors->patchNo)->averageVisibleDirection(scene.GetSurface(surfaceIndex).Normal()).Altitude().degrees() << "\t"
                        << tregenzaSky.GetPatch(factors->patchNo)->averageVisibleDirection(scene.GetSurface(surfaceIndex).Normal()).Azimuth().degrees() << "\t"
                        << theta*180./M_PI << "\t"
                        << factors->unobstructed * groundRadiance << "\t"
                        << ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()->getIAM(theta) << endl;
                        irradiationSW_IAM[surfaceIndex] += factors->unobstructed * groundRadiance * ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()->getIAM(theta);
                    }
#endif
                }
            }
            // diffuse part ground
            irradiationSW[surfaceIndex] += ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_ground() * groundRadiance;
            irradiationDiffuseGround[surfaceIndex] += ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_ground() * groundRadiance;
            // direct part
            if (sunVisibleFarField(pSun->GetPosition().Azimuth().degrees(), pSun->GetPosition().Altitude().degrees()) && Ibn > 0.f) {
                double cosTheta = cosAngleBetween(scene.GetSurface(surfaceIndex).Normal(),pSun->GetPosition());
                if (cosTheta > 0.) {
                    irradiationSW[surfaceIndex] += scene.GetSurface(surfaceIndex).InsolationFactors().GetInsolationFactor(*pSun) * Ibn * cosTheta; // fraction of the surface that is light by the sun
                    irradiationBeam[surfaceIndex] = scene.GetSurface(surfaceIndex).InsolationFactors().GetInsolationFactor(*pSun) * Ibn * cosTheta; // fraction of the surface that is light by the sun
#ifdef DEBUG
                    if (((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()) {
                        float theta = acos(cosTheta);
                        ss_IAM << day << "\t" << hour << "\t" << surfaceIndex << "\t"
                        << scene.GetSurface(surfaceIndex).Normal().Altitude().degrees() << "\t"
                        << scene.GetSurface(surfaceIndex).Normal().Azimuth().degrees() << "\t"
                        << "SUN\t" << "0" << "\t"
                        << pSun->GetPosition().Altitude().degrees() << "\t"
                        << pSun->GetPosition().Azimuth().degrees() << "\t"
                        << theta*180./M_PI << "\t"
                        << scene.GetSurface(surfaceIndex).InsolationFactors().GetInsolationFactor(*pSun) * Ibn * cosTheta << "\t"
                        << ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()->getIAM(theta) << endl;
                        irradiationSW_IAM[surfaceIndex] += scene.GetSurface(surfaceIndex).InsolationFactors().GetInsolationFactor(*pSun) * Ibn * cosTheta
                                                           * ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()->getIAM(theta);
                    }
#endif
                }
            }
        } // end the loop on the surfaces

        // adds the reflections on the surfaces
        irradiationSWn = irradiationSW;
        for (unsigned int r = 0; r < mNbReflections; ++r) {

            // using the CRS sparse matrix format, multiplication of the matrix by the former vector
            vector<float> irradiationSWnew(irradiationSW);
            #pragma omp parallel for schedule(dynamic)
            for (unsigned int i=0; i < getnAi(); ++i) { // loop on the number of elements
                for (unsigned int index=getAi(i); index < getAi(i+1); ++index) {
                    irradiationSWnew[i] += getAn(index) * irradiationSWn[getAj(index)];
                }
            }
            // saves the new irradiation vector in the current
            irradiationSWn = irradiationSWnew;

        }
#ifdef DEBUG
        // adds the inter-reflexions to the irradiance IAM - normal direction efficacy
        float theta = 0.f;
        for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
            if (((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()) {
                irradiationSW_IAM[surfaceIndex] += (irradiationSWn[surfaceIndex] - irradiationSW[surfaceIndex])
                                                    * ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getPVPanel()->getIAM(theta);
            }
        }
#endif

    } // end if sun is up
    else {
        // sets the luminance of the sky and ground to zero
        for (unsigned int patchNo = 0; patchNo < tregenzaSky.PatchCount()/2; ++patchNo) lv[patchNo] = 0.f;
        groundRadiance = 0.f;
        // sets the global irradiance to zero
        pClimate->setIgh(0.f);
    }

    // end of the calculations, save in all surfaces the irradiation
    for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
    {
        // saves all the irradiance
        ((Surface*)(it->SurfaceDelegate()))->setShortWaveIrradiance(irradiationSWn[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

        // saves all the beam irradiance
        ((Surface*)(it->SurfaceDelegate()))->setBeamIrradiance(irradiationBeam[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

        // saves the angle between the sun and the surface
        ((Surface*)(it->SurfaceDelegate()))->setBeamAngle(acosf(GEN::dot_product(it->Normal(),pSun->GetPosition())/(it->Normal().Radius()*pSun->GetPosition().Radius()))/M_PI*180.f);

        // saves all the irradiance IAM
        if (((Surface*)(scene.GetSurface(scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())).SurfaceDelegate()))->getPVPanel()) {
            ((Surface*)(it->SurfaceDelegate()))->setShortWaveIrradiance_IAM(irradiationSW_IAM[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);
#ifdef DEBUG
        ss_IAM_irradiance << day << "\t" << hour << "\t"
                          << scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end()) << "\t"
                          << irradiationSWn[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())] << "\t"
                          << irradiationSW_IAM[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())] << endl;
#endif
        }
    }

    return;

}

void XmlScene::computeDaylight(unsigned int day, unsigned int hour) {

    // initialisation of the sun
    pSun->SetDay(day);
    if (pSun->SetClockTime1(hour) && sky.SetSkyConditions(pClimate->getIdh(day,hour),pClimate->getIbn(day,hour),pSun->GetPosition().Altitude().radians(), pSun->GetPosition().Azimuth().radians(),day)) { // if sun is up and sky conditions are OK

        // start of Igh_vis calculation, computes the global illuminance in the horizontal plane
        float Igh_vis_part = 0.f;
        // direct part contribution to Igh_vis, if sun is above horizon
        if (sunVisibleFarField(pSun->GetPosition().Azimuth().degrees(), pSun->GetPosition().Altitude().degrees())) {
            float cosTheta = GEN::dot_product(GENPoint::Cartesian(0.f,0.f,1.f),pSun->GetPosition());
            if (cosTheta > 0.f) {
                Igh_vis_part += pClimate->getIbn(day,hour) * cosTheta * sky.GetBeamLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour)); // fraction of the surface that is light by the sun
            }
        }
        // diffuse part
        for (unsigned int patchNo = 0; patchNo < tregenzaSky.PatchCount()/2; ++patchNo) // loop on the sky patches
        {
            Igh_vis_part += tregenzaSky.GetPatch(patchNo)->formFactor(GENPoint::Cartesian(0.f,0.f,1.f)) * lv[patchNo] * sky.GetDiffuseLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour));
        }
        pClimate->setIgh_vis(Igh_vis_part);
        // end of Igh_vis calculation

        // illuminance of the surfaces
        vector<float> illuminance(scene.SurfaceCount(), 0.f), illuminanceN(scene.SurfaceCount(), 0.f);

        // use the luminous efficacy to go from irradiance to illuminance
        for (unsigned int i=0; i<scene.SurfaceCount(); ++i)
            illuminance[i] = irradiationBeam[i] * sky.GetBeamLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour)) + (irradiationDiffuseSky[i]+irradiationDiffuseGround[i]) * sky.GetDiffuseLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour));

        // adds the inter-reflections on the surfaces
        illuminanceN = illuminance;
        for (unsigned int r = 0; r < mNbReflections; ++r) {

            // using the CRS sparse matrix format, multiplication of the matrix by the former vector
            vector<float> illuminanceNew(illuminance);
            #pragma omp parallel for schedule(dynamic)
            for (unsigned int i=0; i < getnAi(); ++i) { // loop on the number of elements
                for (unsigned int index=getAi(i); index < getAi(i+1); ++index) {
                    illuminanceNew[i] += getAn(index) * illuminanceN[getAj(index)];
                }
            }
            // saves the new irradiation vector in the current
            illuminanceN = illuminanceNew;

        }

        // end of the calculations, save in all surfaces the illuminance
        for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
        {
            // saves all the irradiation
            ((Surface*)(it->SurfaceDelegate()))->setIlluminance(illuminanceN[scene.SurfaceCount()-distance(it,scene.GetAllSurfaces().end())]);

        }

        // loop on building's surfaces
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
            // only on building surfaces
            if (scene.GetSurface(surfaceIndex).IsBuildingSurface()) {

                // stores the illuminances computed for the two points inside the room
                float internalIlluminance[2];
                // loop on the two measurement points
                for (unsigned int pointIndex=0; pointIndex<2; ++pointIndex) {
                    #ifdef CDEBUG
                    // stores the sky luminances for the actual point
                    float lvSky[tregenzaSky.PatchCount()/2];
                    for (unsigned int patchIndex=0;patchIndex<tregenzaSky.PatchCount()/2;++patchIndex) lvSky[patchIndex] = 0.f;
                    #endif
                    // daylight for point (pointIndex)
                    internalIlluminance[pointIndex] = 0.f;
                    for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).DaylightViewFactors(pointIndex).GetVFs();
                        factors!=scene.GetSurface(surfaceIndex).DaylightViewFactors(pointIndex).GetLastVF();
                        ++factors) // loop on the patches that are non-zero
                    {
                        // for that surface it, output of the informations in the sparse format
                        //logStream << "Patch number: " << factors->patchNo << endl << flush;

                        // the factors->unobstructed contains the patchSolidAngle and the cos of the angles between the surface normal and the center of patch
                        if ( factors->patchNo < static_cast<int>(tregenzaSky.PatchCount()/2) ) {
                            // diffuse part sky and external buildings' luminance
                            internalIlluminance[pointIndex] += factors->unobstructed * lv[factors->patchNo] * sky.GetDiffuseLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour));
                            internalIlluminance[pointIndex] += factors->obstructed * illuminanceN[factors->mainobstructing] * ((Surface*) scene.GetSurface(factors->mainobstructing).SurfaceDelegate())->getShortWaveReflectance() / M_PI;
                            #ifdef CDEBUG
                            lvSky[factors->patchNo] = ( factors->unobstructed * lv[factors->patchNo] * sky.GetDiffuseLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour))
                                                        + factors->obstructed * illuminanceN[factors->mainobstructing] * ((Surface*) scene.GetSurface(factors->mainobstructing).SurfaceDelegate())->getShortWaveReflectance() / M_PI)
                                                         / (factors->unobstructed+factors->obstructed);
                            #endif
                        }
                        // no diffuse part ground
                        // no direct part, perfect imaginary shading device that occludes only the sun
                    }
                    #ifdef CDEBUG
                    // saves the sky luminance in a file (without inter-reflections), for verification only of the shape of the sky
                    exportSkyRadFile(inputFile.substr(0,inputFile.size()-4) + "_SkyLum" + toString(pointIndex) + "_S" + toString(surfaceIndex) + "_" + toString((day-1)*24+hour) + ".rad", lvSky);
                    #endif
                    // computes only a fraction of the incoming flux that goes inside the room
                    internalIlluminance[pointIndex] *= ((Surface*) scene.GetSurface(surfaceIndex).SurfaceDelegate())->getGlazingGvalue();
                }

                // calculation of the internal inter-reflection in the standard room (integrating sphere approximation)
                float downwardIlluminance,upwardIlluminance;
                // adds the diffuse sky component (the direct is considered blocked by solar protections)
                downwardIlluminance = irradiationDiffuseSky[surfaceIndex]*sky.GetDiffuseLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour));
                               //irradiationBeam[surfaceIndex]*sky.GetBeamLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour))
                upwardIlluminance   = irradiationDiffuseGround[surfaceIndex]*sky.GetDiffuseLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour));

                // adds the inter-reflected component
                for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).SWViewFactors().GetVFs();
                    factors!=scene.GetSurface(surfaceIndex).SWViewFactors().GetLastVF();
                    ++factors) // loop on the patches that are non-zero
                {
                    // for that surface it, output of the informations in the sparse format
                    //logStream << "Patch number: " << factors->patchNo << endl << flush;

                    // the factors->unobstructed contains the patchSolidAngle and the cos of the angles between the surface normal and the center of patch
                    if ( factors->patchNo < static_cast<int>(tregenzaSky.PatchCount()/2) ) {
                        // diffuse part sky
                        downwardIlluminance += factors->obstructed * illuminanceN[factors->mainobstructing] * ((Surface*) scene.GetSurface(factors->mainobstructing).SurfaceDelegate())->getShortWaveReflectance() / M_PI;
                    }
                    else {
                        // diffuse part ground
                        upwardIlluminance   += factors->obstructed * illuminanceN[factors->mainobstructing] * ((Surface*) scene.GetSurface(factors->mainobstructing).SurfaceDelegate())->getShortWaveReflectance() / M_PI;
                    }
                }
                // computes the inter-reflected
                float interreflected = 0.f;
                if ( float g = ((Surface*) scene.GetSurface(surfaceIndex).SurfaceDelegate())->getGlazingRatio() > 0.f ) {
                    // computes the reflectances
                    float lowerReflectance = (1.8+5.*(1.-g))/(6.+10.*(1.-g));
                    float upperReflectance = (4.2+3.5*g)/(6.+7.*g);
                    float meanReflectance = 0.5f;
                    float ArOverAw = (22.-3.*g)/3.*g;
                    // computes what is inter-reflected
                    interreflected = ((Surface*) scene.GetSurface(surfaceIndex).SurfaceDelegate())->getGlazingGvalue()*
                                     (downwardIlluminance*lowerReflectance + upwardIlluminance*upperReflectance)/
                                     (1. + ArOverAw*(1.-meanReflectance));
                }

                // saves the illuminance and internally reflected component of the two viewpoints
                ((Wall*) scene.GetSurface(surfaceIndex).SurfaceDelegate())->setInternalIlluminance(internalIlluminance);
                ((Wall*) scene.GetSurface(surfaceIndex).SurfaceDelegate())->setInternallyReflectedIlluminance(interreflected);
            }
        }
    }
    else {
        // no illuminance outside
        pClimate->setIgh_vis(0.f);
        for (DATASurfaceIterator it=scene.GetAllSurfaces();!it.isAtEnd();++it)
            ((Surface*)(it->SurfaceDelegate()))->setIlluminance(0.f);
        // no illuminance in the room
        // loop on building's surfaces
        for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {
            // only on building surfaces
            if (scene.GetSurface(surfaceIndex).IsBuildingSurface()) {
                float internalIlluminance[2];
                internalIlluminance[0] = 0.f;
                internalIlluminance[1] = 0.f;
                #ifdef CDEBUG
                // stores the sky luminances for the actual point
                float lvSky[tregenzaSky.PatchCount()/2];
                for (unsigned int patchIndex=0;patchIndex<tregenzaSky.PatchCount()/2;++patchIndex) lvSky[patchIndex] = 0.f;
                // saves the sky luminance in a file
                for (unsigned int pointIndex=0; pointIndex<2; ++pointIndex)
                    exportSkyRadFile(inputFile.substr(0,inputFile.size()-4) + "_SkyLum" + toString(pointIndex) + "_S" + toString(surfaceIndex) + "_" + toString((day-1)*24+hour) + ".rad", lvSky);
                #endif
                // saves the illuminance of the two viewpoints
                ((Wall*) scene.GetSurface(surfaceIndex).SurfaceDelegate())->setInternalIlluminance(internalIlluminance);
                ((Wall*) scene.GetSurface(surfaceIndex).SurfaceDelegate())->setInternallyReflectedIlluminance(0.f);
            }
        }
    }

    #ifdef CDEBUG
    // saving the luminance in the output file
    exportSkyAndGround(inputFile.substr(0,inputFile.size()-4) + "_SkyRad_" + toString((day-1)*24+hour) + ".rad");
    exportSkyAndGround(inputFile.substr(0,inputFile.size()-4) + "_SkyLum_" + toString((day-1)*24+hour) + ".rad", sky.GetDiffuseLumEffy(pSun->GetPosition().Altitude().radians(), pClimate->getTd(day,hour)));
    #endif

    // end of the method
    return;

}

void XmlScene::computeLongWave(unsigned int day, unsigned int hour) {

    // get the input needed for the longwave balance for the current time-step
    float Ta = pClimate->getToutCelsius(day,hour) + 273.15f; // sky temperature (K)
    float Tg = pClimate->getTgroundCelsius(day,hour) + 273.15f; // ground temperature (K)
    float epsilon_sky = pClimate->getEpsilon_sky(day,hour); // epsilon sky

    // loop on building's surfaces
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int surfaceIndex=0; surfaceIndex<scene.SurfaceCount(); ++surfaceIndex) {

        // the average temperature of the surroundings in kelvins, for this surfaceIndex
        double Tenv4 = 0.;

        // diffuse part sky, use the sum over the visible sky
        Tenv4 += ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_sky()*epsilon_sky*pow(Ta,4);
        // diffuse part ground, use the sum over the visible ground
        Tenv4 += ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_ground()*1.0*pow(Tg,4); // emissivity of the ground is close to 1 according to Meteonorm Reference

        // loop on the Tregenza patches seen from surface(surfaceIndex)
        for (DATAViewFactorSetSparse::const_iterator factors=scene.GetSurface(surfaceIndex).SWViewFactors().GetVFs();
            factors!=scene.GetSurface(surfaceIndex).SWViewFactors().GetLastVF();
            ++factors)
        {
            // obstructed part of the Tregenza patch -> meaning obstructed by a surface (factors->mainobstructing)
            if (factors->obstructed > 0.f) {

                //logStream << "MainObstructing: " << factors->mainobstructing << " has a temperature defined." << endl << flush;
                float Tmainobstructing = ((Surface*)(scene.GetSurface(factors->mainobstructing).SurfaceDelegate()))->getTemperature() + 273.15f;
                float emissivity_mainobstructing = ((Surface*)(scene.GetSurface(factors->mainobstructing).SurfaceDelegate()))->getLongWaveEmissivity();

                if (!isnan(Tmainobstructing)) {
                    //cout << "Tmainobstructing=" << Tmainobstructing << endl;
                    // check physical limits
                    if (Tmainobstructing > 500.f) throw string("Main obstructing surface id(key): " + toString(((Surface*)(scene.GetSurface(factors->mainobstructing).SurfaceDelegate()))->getId())
                                                               + "(" + toString(((Surface*)(scene.GetSurface(factors->mainobstructing).SurfaceDelegate()))->getKey()) + ")"
                                                               +", Temperature: " + toString(Tmainobstructing));
                    // compute the average temperature^4
                    Tenv4 += factors->obstructed*emissivity_mainobstructing*pow(Tmainobstructing,4);
                }
                else {/* we make the assumption that the obtructing surface is another building whose temperature is unknown
                         and supposed to be the same as the one from this surfaceIndex
                         -> we don't take any temperature and don't take into account that factors->obstructed */
                }
            }
        }

        // we divide by the projected solid angle of the hemisphere as we want an equivalent radiating temperature for the whole hemisphere
        Tenv4 /= ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_hemisphere();
        //logStream << "PI approx. " << ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->getProjectedSolidAngle_hemisphere() << endl << flush;

        // no loop for inter-reflections here, only average temperature used
        // saves the environmental temperature un the Surface itself
        ((Surface*)(scene.GetSurface(surfaceIndex).SurfaceDelegate()))->setEnvironmentalTemperature(pow(Tenv4,0.25)-273.15);
    }

    return;

}

int XmlScene::computeWarmUp() {

    vector<int> value(pDistrict->getnBuildings(),0);
    #pragma omp parallel for schedule(dynamic)
    for (size_t i=0; i<pDistrict->getnBuildings(); ++i){
        value[i] = Model::ThermalWarmUpTime(pDistrict->getBuilding(i));
    }
    if (pDistrict->getnBuildings()==0) return 0;
    else return *max_element(value.begin(),value.end());

}

void XmlScene::computeThermal(unsigned int day, unsigned int hour) {

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int i=0; i<pDistrict->getnBuildings(); ++i) {
        //logStream << "Building " << i << "\tNumber of thread: " << omp_get_thread_num() << endl << flush;
        Model::ThermalStepImplicit(pDistrict->getBuilding(i),pClimate,day,hour);
    }

    // case of an energy hub connected to many buildings
    Model::noHVAC_Control_EnergyHub(pDistrict,pClimate,day,hour);


    // evaluates the new buildings' temperature according to what is available
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int i=0; i<pDistrict->getnBuildings(); ++i) {
        if (Model::thermalExplicit) {
            Model::ThermalStepExplicitTemperature(pDistrict->getBuilding(i),pClimate,day,hour);
            timeSteps2Simulated+=Model::dt/Model::dt2;
        }
        else Model::ThermalStepImplicitTemperature(pDistrict->getBuilding(i),pClimate,day,hour);
    }

    // evaluates the tree temperature
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int i=0;i<pDistrict->getnTrees();++i) {
        Model::ThermalStepTree(pDistrict->getTree(i),pClimate,day,hour);
    }

    // evaluates the ground temperature
    #pragma omp parallel
    {
        for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
            #pragma omp single nowait
            {
                if ((*it)->isDetailedSimulation())
                    Model::ThermalStepImplicitTemperature((*it),pClimate,day,hour);
                else
                    Model::ThermalStepImplicitTemperature_simplified((*it),pClimate,day,hour);
            }
        }
    }

    return;
}

#ifdef FMI
void XmlScene::initialiseThermal_EnergyPlus() {

    // initialise the FMU
    fmi1_callback_functions_t callBackFunctions;
    const char* FMUPath;
    const char* tmpPath;
    jm_callbacks callbacks;
    fmi_version_enu_t version;

    // loop on all the buildings simulated with FMU
    for (unsigned int i=0; i<pDistrict->getnBuildings(); ++i) {
        // test of the building is simulated with E+
        if (pDistrict->getBuilding(i)->isEP()) {

            // gets the paths and name of the EnergyPlus FMU file
            FMUPath = pDistrict->getBuilding(i)->getFMUFile().c_str();
            tmpPath = pDistrict->getBuilding(i)->getTMPPath().c_str();

            callbacks.malloc = malloc;
            callbacks.calloc = calloc;
            callbacks.realloc = realloc;
            callbacks.free = free;
            callbacks.logger = importlogger;
            callbacks.log_level = jm_log_level_debug;
            callbacks.context = 0;

            callBackFunctions.logger = fmi1_log_forwarding;
            callBackFunctions.allocateMemory = calloc;
            callBackFunctions.freeMemory = free;

            context = fmi_import_allocate_context(&callbacks);

            // this actually unzipts the fmu to tmpPath and then reads in the xml file
            version = fmi_import_get_fmi_version(context, FMUPath, tmpPath);

            if (version != fmi_version_1_enu)
            {
                throw(string("Only version 1.0 is supported so far\n"));
            }

            fmu = fmi1_import_parse_xml(context, tmpPath);

            if (!fmu)
            {
                throw(string("Error parsing XML, exiting\n"));
            }

            jmstatus = fmi1_import_create_dllfmu(fmu, callBackFunctions, 1); // FIXME: try setting registerGlobally to 0 for thread safety
            if (jmstatus == jm_status_error)
            {
                throw(string("Could not create the DLL loading mechanism(C-API) (error:" + string(fmi1_import_get_last_error(fmu)) + ").\n"));
            }

            // initializing the EnergyPlus simulation
            fmi1_string_t instanceName = "RevitToCitySim";
            fmi1_string_t fmuLocation = 0;
            fmi1_string_t mimeType = "";
            fmi1_real_t timeout = 0.0;
            fmi1_boolean_t visible = fmi1_false;
            fmi1_boolean_t interactive = fmi1_false;

            // Simulation setup: Simulation starts at midnight (this cannot be changed)
            fmi1_real_t tstart = (beginDay-1)* 24.0 * 60 * 60; // starts usually at 0.0
            fmi1_real_t tend = endDay * 24.0 * 60 * 60; // 365 days
            fmi1_boolean_t StopTimeDefined = fmi1_false;

            jmstatus = fmi1_import_instantiate_slave(fmu, instanceName, fmuLocation, mimeType, timeout, visible, interactive);
            if (jmstatus == jm_status_error)
            {
                throw(string("fmi1_import_instantiate_slave failed\n"));
            }

            // this line starts up EnergyPlus and does the warm up
            fmistatus = fmi1_import_initialize_slave(fmu, tstart, StopTimeDefined, tend);
            if (fmistatus != fmi1_status_ok)
            {
                throw(string("fmi1_import_initialize_slave failed\n"));
            }

        }

    } // end of the loop on Buildings

}

fmi1_value_reference_t XmlScene::getValueReference(string variableName) {

    fmi1_import_variable_t* pVariable = fmi1_import_get_variable_by_name(fmu,variableName.c_str());
    if (pVariable == NULL) { throw(string("Error: Variable " + variableName + " not found in EnergyPlus simulation.")); }
    return fmi1_import_get_variable_vr(pVariable);

}

void XmlScene::computeThermal_EnergyPlus(unsigned int day, unsigned int hour) {

    // evaluates the new buildings' temperature according to what is available
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int i=0; i<pDistrict->getnBuildings(); ++i) {
        if (pDistrict->getBuilding(i)->isEP()) {

            // defines the values to take out of the simulation
            fmi1_real_t realValue;
            fmi1_value_reference_t valueReference;

            // the buffer of string
            string buffer;

            // Note: you can also set values here! But they need to be set up in the IDF as well...
            buffer = "Outdoor Drybulb";
                    valueReference = getValueReference(buffer);
            realValue = pClimate->getToutCelsius(day,hour);
            fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);

            buffer = "Outdoor Dewpoint";
                    valueReference = getValueReference(buffer);
            realValue = pClimate->getTd(day,hour);
            fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);

            buffer = "Outdoor Relative Humidity";
                    valueReference = getValueReference(buffer);
            realValue = pClimate->getRelativeHumidity(day,hour);
            fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);

            buffer = "Diffuse Solar";
                    valueReference = getValueReference(buffer);
            realValue = pClimate->getIdh(day,hour);
            fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);

            buffer = "Direct Solar";
                    valueReference = getValueReference(buffer);
            realValue = pClimate->getIbn(day,hour);
            fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);

            buffer = "Wind Speed";
                    valueReference = getValueReference(buffer);
            realValue = pClimate->getWindSpeed(day,hour);
            fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);

            buffer = "Wind Direction";
                    valueReference = getValueReference(buffer);
            realValue = pClimate->getWindDirection(day,hour);
            fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);

            for (unsigned int j=0; j<pDistrict->getBuilding(i)->getnZones(); ++j) {
                buffer = pDistrict->getBuilding(i)->getZone(j)->getEp_id() + "::Occupation";
                    valueReference = getValueReference(buffer);
                realValue = pDistrict->getBuilding(i)->getZone(j)->getOccupantsCount();
                fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);
            }

            // sets the variable for the LW exchange
            for (unsigned int j=0; j<pDistrict->getBuilding(i)->getnZones(); ++j) {
                // set the last time step results in the walls
                for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnWalls(); ++k) {
                    // sets Tenv
                    buffer = pDistrict->getBuilding(i)->getZone(j)->getWall(k)->getEp_id() + "::Tenv";
                    valueReference = getValueReference(buffer);
                    realValue = pDistrict->getBuilding(i)->getZone(j)->getWall(k)->getEnvironmentalTemperature();
                    fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);
                    // sets Hr
                    buffer = pDistrict->getBuilding(i)->getZone(j)->getWall(k)->getEp_id() + "::Henv";
                    valueReference = getValueReference(buffer);
                    realValue = pDistrict->getBuilding(i)->getZone(j)->getWall(k)->get_hr();
                    fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);
                }
                // sset the last time step results in the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnRoofs(); ++k) {
                    // sets Tenv
                    buffer = pDistrict->getBuilding(i)->getZone(j)->getRoof(k)->getEp_id() + "::Tenv";
                    valueReference = getValueReference(buffer);
                    realValue = pDistrict->getBuilding(i)->getZone(j)->getRoof(k)->getEnvironmentalTemperature();
                    fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);
                    // sets Hr
                    buffer = pDistrict->getBuilding(i)->getZone(j)->getRoof(k)->getEp_id() + "::Henv";
                    valueReference = getValueReference(buffer);
                    realValue = pDistrict->getBuilding(i)->getZone(j)->getRoof(k)->get_hr();
                    fmistatus = fmi1_import_set_real(fmu, &valueReference, 1, &realValue);
                }
            }
            // end of the loop on Zones

            // simulate with EnergyPlus
            fmi1_boolean_t newStep = fmi1_true;
            fmi1_real_t tcur = (((day-1)%365 * 24) + hour - 1) * 60 * 60; // remain within the boundaries 1-365, JK - 21.06.2015
            fmi1_real_t hstep = 60 * 60; // one hour, must match timestep in IDF, so for TIMESTEP = 4, use 15*60

            // call fmi1_import_do_step for each TIMESTEP in the simulation
            fmistatus = fmi1_import_do_step(fmu, tcur, hstep, newStep);

            // saves the result in the zones
            for (unsigned int j=0; j<pDistrict->getBuilding(i)->getnZones(); ++j) {
                // gets values from the E+ simulation
                buffer = pDistrict->getBuilding(i)->getZone(j)->getEp_id() + "::Total Heating Energy";
                    valueReference = getValueReference(buffer);
                fmistatus = fmi1_import_get_real(fmu, &valueReference, 1, &realValue);
                // sets the heating needs
                pDistrict->getBuilding(i)->getZone(j)->eraseHeating_back();
                pDistrict->getBuilding(i)->getZone(j)->setHeating(realValue/3600.);
                // gets values from the E+ simulation
                buffer = pDistrict->getBuilding(i)->getZone(j)->getEp_id() + "::Total Cooling Energy";
                    valueReference = getValueReference(buffer);
                fmistatus = fmi1_import_get_real(fmu, &valueReference, 1, &realValue);
                // sets the cooling needs
                pDistrict->getBuilding(i)->getZone(j)->eraseCooling_back();
                pDistrict->getBuilding(i)->getZone(j)->setCooling(-realValue/3600.);
                // gets the ventilation flow rate
                buffer = pDistrict->getBuilding(i)->getZone(j)->getEp_id() + "::Ventilation Volume Flow Rate";
                    valueReference = getValueReference(buffer);
                fmistatus = fmi1_import_get_real(fmu, &valueReference, 1, &realValue);
                // sets the ventilation rate
                pDistrict->getBuilding(i)->getZone(j)->eraseVdotVent_back();
                pDistrict->getBuilding(i)->getZone(j)->setVdotVent(realValue);
                // gets the ambient air temperature
                buffer = pDistrict->getBuilding(i)->getZone(j)->getEp_id() + "::Mean Air Temperature";
                    valueReference = getValueReference(buffer);
                fmistatus = fmi1_import_get_real(fmu, &valueReference, 1, &realValue);
                // saves the ambient air temperature
                pDistrict->getBuilding(i)->getZone(j)->eraseT_back();
                pDistrict->getBuilding(i)->getZone(j)->setTa(realValue);
                // saves the results in the walls
                for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnWalls(); ++k) {
                    // gets values from the E+ simulation
                    buffer = pDistrict->getBuilding(i)->getZone(j)->getWall(k)->getEp_id() + "::Outside Surface Temperature";
                    valueReference = getValueReference(buffer);
                    fmistatus = fmi1_import_get_real(fmu, &valueReference, 1, &realValue);
                    // sets the values
                    pDistrict->getBuilding(i)->getZone(j)->getWall(k)->eraseTemperature_back();
                    pDistrict->getBuilding(i)->getZone(j)->getWall(k)->setTemperature(realValue);
                }
                // saves the results in the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(i)->getZone(j)->getnRoofs(); ++k) {
                    // gets values from the E+ simulation
                    buffer = pDistrict->getBuilding(i)->getZone(j)->getRoof(k)->getEp_id() + "::Outside Surface Temperature";
                    valueReference = getValueReference(buffer);
                    fmistatus = fmi1_import_get_real(fmu, &valueReference, 1, &realValue);
                    // sets the values
                    pDistrict->getBuilding(i)->getZone(j)->getRoof(k)->eraseTemperature_back();
                    pDistrict->getBuilding(i)->getZone(j)->getRoof(k)->setTemperature(realValue);
                }
            // end of the loop on Zones
            }
        }
    }

    return;
}

void XmlScene::terminateThermal_EnergyPlus() {

    for (unsigned int i=0; i<pDistrict->getnBuildings(); ++i) {
        if (pDistrict->getBuilding(i)->isEP()) {

            fmistatus = fmi1_import_terminate_slave(fmu);
            if (fmistatus != fmi1_status_ok)
            {
                throw(string("fmi1_import_terminate_slave failed\n"));
            }
            fmi1_import_free_slave_instance(fmu);
            // clean up
        //    fmi1_import_destroy_dllfmu(fmu);
        //    fmi1_import_free(fmu);
        //    fmi_import_free_context(context);

            return;
        }
    }

}
#endif

void XmlScene::simulate() {

    // initialises the Far Field obstructions vector
    initialiseFarField();

    // gives the number of threads available
    #if !defined(DEBUG) && !defined(MONOTHREAD)
    logStream << "OpenMP max threads: " << omp_get_max_threads() << endl << flush;
    #endif
    //#ifdef WS
    //omp_set_num_threads(min(4,omp_get_max_threads()));
    //#endif

    // initialises the FMU
    #ifdef FMI
    initialiseThermal_EnergyPlus();
    #endif

    initialiseSimulationParameters();

    updateSimulationModelParameters();

    // pre-conditioning period
    int warmUpDays = computeWarmUp();
    logStream << "Pre-conditioning period: " << warmUpDays << " days" << endl;
    for (int day = -warmUpDays; day <= -1; ++day) {

        unsigned int preDay = (day+beginDay+364)%365+1;
        logStream << "Pre-conditioning period - day: " << preDay << endl << flush;
        for (unsigned int hour = 1; hour <= 24; ++hour) {

            simulateTimeStep(preDay,hour,true);
        }
    }

    #ifdef DEBUG
    writeTHResultsText(getInputFileNoExt() + "_TH.out");
    #endif // DEBUG

    if(warmUpDays>0) clearPreConditioningResults();

    // debut de la simulation sur la periode consideree
    for (unsigned int day = beginDay; day<=endDay; ++day) {

        logStream << "Day: " << day << "\tMemory usage: " << setprecision(4) << static_cast<float>(memoryUsage())/1.e6f << " Mb" << endl << flush;
        for (unsigned int hour = 1; hour <= 24; ++hour) {

            simulateTimeStep(day,hour,false);
        }

        checkMemoryUsage(0.5e9 ); // save intermediate results and clear if greater than 500 Mb
    }

    logStream << "Simulation ended." << endl << flush;
    #ifdef FMI
    terminateThermal_EnergyPlus();
    #endif

}

void XmlScene::simulateCumulativeIrradiance() {

    // initialises the Far Field obstructions vector
    initialiseFarField();

    // gives the number of threads available
    #if !defined(DEBUG) && !defined(MONOTHREAD)
    logStream << "OpenMP max threads: " << omp_get_max_threads() << endl << flush;
    #endif

    initialiseSimulationParameters();

    // computes the annual cumulative diffuse sky
    computeCumulativeRadiance(beginDay,endDay);
    computeFarField();

    // computes the annual cumulative irradiation
    computeShortWave_Diffuse();

    // computes the sum of the hourly beam irradiance
    for (unsigned int day = beginDay; day<=endDay; day++) {
        logStream << "Day: " << day << endl;
        for (unsigned int hour = 1; hour <= 24; hour++) {
            computeShortWave_Beam(day,hour);
        }
    }

    // computes the interreflected component
    computeShortWave_Interreflected();

    // increase the number of simulated time steps
    ++timeStepsSimulated;
    logStream << "Simulation ended." << endl;

}

void XmlScene::simulateRadiation(bool daylight) {

    // initialises the Far Field obstructions vector
    initialiseFarField();

    // gives the number of threads available
    #if !defined(DEBUG) && !defined(MONOTHREAD)
    logStream << "OpenMP max threads: " << omp_get_max_threads() << endl << flush;
    #endif

    initialiseSimulationParameters();

    // debut de la simulation sur la periode consideree
    for (unsigned int day = beginDay; day<=endDay; day++) {

        logStream << "Day: " << day << "\tMemory usage: " << setprecision(4) << static_cast<float>(memoryUsageIrradiation())/1.e6f << " Mb" << endl << flush;
        for (unsigned int hour = 1; hour <= 24; hour++) {

            simulateTimeStep(day, hour, false, true, daylight);
        }

        checkMemoryUsage(0.5e9 ,/*radOnly=*/ true);// save intermediate results and clear if greater than 500 Mb
    }
    logStream << "Simulation ended." << endl << flush;

}

void XmlScene::updateSimulationModelParameters(){

    // ground surfaces
    for (forward_list<Ground*>::iterator itSurfaces = getDistrict()->getGrounds()->begin(); itSurfaces != getDistrict()->getGrounds()->end(); ++itSurfaces) {
        (*itSurfaces)->initialiseModel();
    }

    // Buildings
    for (vector<Building*>::iterator itBuildings = getDistrict()->getBuildings()->begin();
         itBuildings != getDistrict()->getBuildings()->end();
         ++itBuildings) (*itBuildings)->update();

}

void XmlScene::initialiseSimulationParameters() {

    eraseResults(0,true);
    timeStepsSimulated=0;
    timeSteps2Simulated=0;
    simulationIndex=0;
    preTimeStepsSimulated=0;
    preTimeSteps2Simulated=0;

}

void XmlScene::simulateTimeStep(int day, int hour, bool preCond, bool radOnly, bool doDayLightSim){
    // output of the current moment
    // logStream << "Day: " << day << "\tHour: " << hour << endl << flush;
    // logStream << "computeShortWave" << endl;
    computeShortWave(day,hour);

    // exportSWFile( inputFile.substr(0,inputFile.size()-4) + "_SW_" + toString(day) + "_" + toString(hour) + ".out" );
    // logStream << "computeDaylight" << endl;
    if (doDayLightSim) computeDaylight(day,hour);

    if (!radOnly) {
        //logStream << "computeThermal" << endl;
        computeThermal(day,hour);
        if (!preCond){
            #ifdef FMI
            computeThermal_EnergyPlus(day,hour);
            #endif
        }
        //logStream << "computeLongWave" << endl;
        computeLongWave(day,hour);
    }

    if (preCond){
        ++preTimeStepsSimulated;
    }
    else {
        ++timeStepsSimulated;
    }
}

void XmlScene::clearPreConditioningResults(){
    // erasing of the results provided by the pre-conditioning period, keeping only one value (last time step)
    eraseResults(1,true);
    preTimeStepsSimulated = 1;
    if (Model::thermalExplicit) preTimeSteps2Simulated = 1;
}

void XmlScene::checkMemoryUsage(double limit, bool radOnly){
    // check the memory content
    if ( radOnly && memoryUsage() > limit ) { // 0.5e9 -> greater than 500 Mb
        logStream << "Saving results..." << endl << flush;
        writeSWResultsText(inputFile.substr(0,inputFile.size()-4) + "_SW.out");
        //writeDLResultsText(inputFile.substr(0,inputFile.size()-4) + "_DL.out");
        // zeroing the counters, and saving the current index of the simulation
        eraseResultsIrradiation(1);
        simulationIndex += timeStepsSimulated;
        timeStepsSimulated = 0;
    }
    else if ( memoryUsage() > limit ) { // 0.5e9 -> greater than 500 Mb
        logStream << "Saving results..." << endl << flush;
        #ifndef WS
        writeSWResultsText(inputFile.substr(0,inputFile.size()-4) + "_SW.out");
        writeSWvResultsText(inputFile.substr(0,inputFile.size()-4) + "_SWv.out");
        writeDLResultsText(inputFile.substr(0,inputFile.size()-4) + "_DL.out");
        writeLWResultsText(inputFile.substr(0,inputFile.size()-4) + "_LW.out");
        writeTHResultsText(inputFile.substr(0,inputFile.size()-4) + "_TH.out");
        writeTSResultsText(inputFile.substr(0,inputFile.size()-4) + "_TS.out");
        writeHCResultsText(inputFile.substr(0,inputFile.size()-4) + "_HC.out");
        writeETResultsText(inputFile.substr(0,inputFile.size()-4) + "_ET.out");
        writeCMResultsText(inputFile.substr(0,inputFile.size()-4) + "_CM.out");
        if (Model::thermalExplicit)
            writeTHExplicitResultsText(inputFile.substr(0,inputFile.size()-4) + "_THExplicit.out");
        #endif
        // zeroing the counters, and saving the current index of the simulation
        eraseResults(1,false);
        simulationIndex += timeStepsSimulated;
        timeStepsSimulated = 0;
        timeSteps2Simulated = 0;
    }
}

void XmlScene::eraseResultsThermal_back() {

    // deleting the last result for the TH
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            pDistrict->getBuilding(j)->getZone(k)->eraseT_back();
            // the two following elements do not erase the results if not in preSimulation
            pDistrict->getBuilding(j)->getZone(k)->eraseHeating_back();
            pDistrict->getBuilding(j)->getZone(k)->eraseCooling_back();
            if (pDistrict->getBuilding(j)->getHVACpresence()) {
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACHeat_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACHeatAvailable_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACCool_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACCoolAvailable_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACHumidification_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACHumidificationAvailable_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACEvaporation_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACEvaporationAvailable_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACReheat_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACReheatAvailable_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACMassFlowRate_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseHVACMassFlowRateAvailable_back();
            }
            pDistrict->getBuilding(j)->getZone(k)->eraseQi_back();
            pDistrict->getBuilding(j)->getZone(k)->eraseQs_back();
            pDistrict->getBuilding(j)->getZone(k)->eraseVdotVent_back();
        }
        if (pDistrict->getBuilding(j)->getHeatStock() != NULL)
            pDistrict->getBuilding(j)->eraseHeatStockTemperature_back();
        if (pDistrict->getBuilding(j)->getColdStock() != NULL)
            pDistrict->getBuilding(j)->eraseColdStockTemperature_back();
        pDistrict->getBuilding(j)->eraseMachinePower_back();
        pDistrict->getBuilding(j)->eraseFuelConsumption_back();
        pDistrict->getBuilding(j)->eraseElectricConsumption_back();
        pDistrict->getBuilding(j)->eraseSolarPVProduction_back();
        pDistrict->getBuilding(j)->eraseSolarThermalProduction_back();
    }

//beginning of contents added by Dapeng // Cognet: Adapted to new code.
    for(unsigned int i=0; i<pDistrict->getnDECs(); ++i) {
        pDistrict->getDEC(i)->eraseRecords_back();
//        pDistrict->getDEC(i)->eraseThermalLosses_back();
//        pDistrict->getDEC(i)->erasePumpPowers_back();
//        pDistrict->getDEC(i)->getNetwork()->eraseTemperatures_back();
//        pDistrict->getDEC(i)->getNetwork()->erasePressureDiffs_back();
//        pDistrict->getDEC(i)->getNetwork()->eraseMassFlows_back();
//        pDistrict->getDEC(i)->eraseFuelConsumptions_back();
//        pDistrict->getDEC(i)->eraseElectricConsumptions_back();
    }
//end of contents added by Dapeng

    // deleting the results for the THExplicit
    // emptying the vectors (only leaving one value per vector)
    if (Model::thermalExplicit) {
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            // loop on the zones
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                pDistrict->getBuilding(j)->getZone(k)->eraseTExpl_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseWindowState_back();
                pDistrict->getBuilding(j)->getZone(k)->eraseVdotVent_back();
                for (unsigned int l=0; l<pDistrict->getBuilding(j)->getZone(k)->getnWalls(); ++l) {
                    pDistrict->getBuilding(j)->getZone(k)->getWall(l)->eraseLowerShadingState_back();
                    pDistrict->getBuilding(j)->getZone(k)->getWall(l)->eraseUpperShadingState_back();
                }
            }
        }
    }

    // deleting the results for the TS
    // emptying the vectors (leaving only one value in each vector)
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseTemperature_back();
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseTemperature_back();
            }
            // loop on the obstructing surfaces
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseTemperature_back();
            }
        } // end the loop on zones
    }
    // loop on the tree surfaces in the district itself
    for (size_t j=0; j<pDistrict->getnTrees(); ++j)
        for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
            pDistrict->getTree(j)->getSurface(k)->eraseTemperature_back();

    // emptying the vectors of the ground temperature
    for (forward_list<Ground*>::iterator it=getDistrict()->getGrounds()->begin();it!=getDistrict()->getGrounds()->end();++it) {
        (*it)->eraseTemperature_back();
    }

    return;

}

/**
 * @brief XmlScene::eraseResults
 * @param keepValue : The number of values to be left in the vectors (usually 0 before simulation, 1 after pre-process)
 * @param eraseAllResults : If false erase only intermediate results, leaving main results available
 */
void XmlScene::eraseResults(unsigned int keepValue, bool eraseAllResults) {

    //logStream << "erase results " << endl;

    // deleting the results for the global horizontal
    // emptying the vectors (leaving only keepValue in each vector)
    if (pClimate!=NULL) {
        pClimate->eraseIgh(keepValue);
        pClimate->eraseIgh_vis(keepValue);
    }

    // deleting all results, loop on buildings
    // emptying the vectors (leaving only one value in each vector)
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        //logStream << "Loop buildings: " << pDistrict->getBuilding(j)->getKey() << endl;

        // deleting the results for the TH
        if (pDistrict->getBuilding(j)->getHeatStock() != NULL)
            pDistrict->getBuilding(j)->eraseHeatStockTemperature(keepValue);
        if (pDistrict->getBuilding(j)->getDHWHeatStock() != NULL)
            pDistrict->getBuilding(j)->eraseDHWStockT(keepValue);
        if (pDistrict->getBuilding(j)->getColdStock() != NULL)
            pDistrict->getBuilding(j)->eraseColdStockTemperature(keepValue);
        if (eraseAllResults){
            pDistrict->getBuilding(j)->eraseMachinePower(keepValue); // DP: modified to keep MachinePower results (MEU webservice)
            pDistrict->getBuilding(j)->eraseFuelConsumption();     // DP: Probably not necessary for MEU
            pDistrict->getBuilding(j)->eraseElectricConsumption(); // DP: Probably not necessary for MEU
            pDistrict->getBuilding(j)->eraseSolarPVProduction();
            pDistrict->getBuilding(j)->eraseSolarThermalProduction();
        }

        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {

            // deleting the results for the TH
            pDistrict->getBuilding(j)->getZone(zone)->eraseT(keepValue);
            // the two following elements do not erase the results if not in preSimulation
            if(eraseAllResults){ // && keepValue???
                pDistrict->getBuilding(j)->getZone(zone)->eraseHeating();
                pDistrict->getBuilding(j)->getZone(zone)->eraseCooling();
            }
            if (pDistrict->getBuilding(j)->getHVACpresence()) {
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACHeat(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACHeatAvailable(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACCool(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACCoolAvailable(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACHumidification(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACHumidificationAvailable(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACEvaporation(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACEvaporationAvailable(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACReheat(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACReheatAvailable(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACMassFlowRate(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->eraseHVACMassFlowRateAvailable(keepValue);
            }
            pDistrict->getBuilding(j)->getZone(zone)->eraseQi(keepValue);
            pDistrict->getBuilding(j)->getZone(zone)->eraseQs(eraseAllResults); // never keep a value in the vector
            pDistrict->getBuilding(j)->getZone(zone)->eraseVdotVent(keepValue);


            // deleting the results for the SW, LW, TS, DL and DLout
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseShortWaveIrradiance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseIlluminance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseEnvironmentalTemperature(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseInternalIlluminance0(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseInternalIlluminance1(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseInternallyReflectedIlluminance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseTemperature(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->erase_hc(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->erasePVElectricProduction(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseSolarThermalProduction(keepValue);
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseShortWaveIrradiance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseIlluminance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseEnvironmentalTemperature(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseTemperature(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->erase_hc(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseWaterEvapotranspiration(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->erasePVElectricProduction(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseSolarThermalProduction(keepValue);
            }
            // loop on the obstructing surfaces
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseShortWaveIrradiance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseIlluminance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseEnvironmentalTemperature(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseTemperature(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->erasePVElectricProduction(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseSolarThermalProduction(keepValue);
            }
        } // end the loop on zones
    }


    // deleting the results for the THExplicit
    // emptying the vectors (only leaving one value per vector)
    if (Model::thermalExplicit) {
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            // loop on the zones
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                pDistrict->getBuilding(j)->getZone(k)->eraseTExpl(keepValue);
                pDistrict->getBuilding(j)->getZone(k)->eraseWindowState(keepValue);
                pDistrict->getBuilding(j)->getZone(k)->eraseVdotVent(keepValue);
                for (unsigned int l=0; l<pDistrict->getBuilding(j)->getZone(k)->getnWalls(); ++l) {
                    pDistrict->getBuilding(j)->getZone(k)->getWall(l)->eraseLowerShadingState(keepValue);
                    pDistrict->getBuilding(j)->getZone(k)->getWall(l)->eraseUpperShadingState(keepValue);
                }
            }
        }
    }
    // loop on the tree surfaces in the district itself
    for (size_t j=0; j<pDistrict->getnTrees(); ++j) {
        for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k) {
            pDistrict->getTree(j)->getSurface(k)->eraseTemperature(keepValue);
            pDistrict->getTree(j)->getSurface(k)->eraseShortWaveIrradiance(keepValue);
            pDistrict->getTree(j)->getSurface(k)->eraseEnvironmentalTemperature(keepValue);
        }
    }
    // emptying the vectors of the surfaces irradiation & temp
    for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j) {
        pDistrict->getSurface(j)->eraseShortWaveIrradiance(keepValue);
        pDistrict->getSurface(j)->eraseEnvironmentalTemperature(keepValue);
    }
    // emptying the vectors of the ground irradiation & temp
    for (forward_list<Ground*>::iterator it=getDistrict()->getGrounds()->begin();it!=getDistrict()->getGrounds()->end();++it) {
        (*it)->eraseTemperature(keepValue);
        (*it)->eraseWaterEvapotranspiration(keepValue);
        (*it)->eraseShortWaveIrradiance(keepValue);
        (*it)->eraseEnvironmentalTemperature(keepValue);
        (*it)->erase_hc(keepValue);
    }

//beginning of contents added by Dapeng // Cognet: Adapted to new code.
    for(unsigned int i=0; i<pDistrict->getnDECs(); ++i) {
        pDistrict->getDEC(i)->eraseRecords(keepValue);
//        pDistrict->getDEC(i)->eraseMachinePowers(keepValue);
//        pDistrict->getDEC(i)->eraseThermalLosses(keepValue);
//        pDistrict->getDEC(i)->erasePumpPowers(keepValue); makeThisSubFunctions
//        pDistrict->getDEC(i)->getNetwork()->eraseTemperatures(keepValue);
//        pDistrict->getDEC(i)->getNetwork()->erasePressureDiffs(keepValue);
//        pDistrict->getDEC(i)->getNetwork()->eraseMassFlows(keepValue);
//        pDistrict->getDEC(i)->eraseFuelConsumptions(keepValue);
//        pDistrict->getDEC(i)->eraseElectricConsumptions(keepValue);
    }
//end of contents added by Dapeng

    return;

}

void XmlScene::eraseResultsIrradiation(unsigned int keepValue) {

    // deleting the results for the SW
    // emptying the vectors (leaving only one value in each vector)
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseShortWaveIrradiance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseIlluminance(keepValue);
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseShortWaveIrradiance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->eraseIlluminance(keepValue);
            }
            // loop on the obstructing surfaces
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseShortWaveIrradiance(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->eraseIlluminance(keepValue);
            }
        } // end the loop on zones
    }

    // loop on the tree surfaces in the district itself
    for (size_t j=0; j<pDistrict->getnTrees(); ++j)
        for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
            pDistrict->getTree(j)->getSurface(k)->eraseShortWaveIrradiance(keepValue);

    // emptying the vectors of the ground irradiation
    for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j) {
        pDistrict->getSurface(j)->eraseShortWaveIrradiance(keepValue);
    }
    // emptying the vectors of the ground irradiation
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        (*it)->eraseShortWaveIrradiance(keepValue);
    }

    // deleting the results for the DL
    // emptying the vectors (leaving only one value in each vector)
    pClimate->eraseIgh(keepValue);
    pClimate->eraseIgh_vis(keepValue);
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseInternalIlluminance0(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseInternalIlluminance1(keepValue);
                pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->eraseInternallyReflectedIlluminance(keepValue);
            }
        } // end the loop on zones
    }

    return;

}

unsigned int XmlScene::getColumnIndex(Surface *surface) {

    unsigned int counter = 1; // the first column is the Igh
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                if (surface == pDistrict->getBuilding(j)->getZone(zone)->getWall(k)) return counter++;
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                if (surface == pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)) return counter++;
            }
        }
    }
    throw(string("Surface id=" + toString(surface->getId()) + " not found."));
    return 0;
}

void XmlScene::writeSWHeaderText(string fileOut, string unit) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#" << "Igh" << "\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" <<
                            pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << "(" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getKey() << "):" << unit << "\t";
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" <<
                            pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getId() << "(" << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getKey() << "):" << unit << "\t";
            }
            // loop on the obstructing surfaces
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" <<
                            pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getId() << "(" << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getKey() << "):" << unit << "\t";
            }
        } // end the loop on zones
    }

    // loop on the tree surfaces in the district itself
    for (size_t j=0; j<pDistrict->getnTrees(); ++j)
        for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
            textFile << pDistrict->getTree(j)->getId() << "(" << pDistrict->getTree(j)->getKey() << "):" <<
                        pDistrict->getTree(j)->getSurface(k)->getId() << "(" << pDistrict->getTree(j)->getSurface(k)->getKey() << "):" << unit << "\t";

    // loop on the obstructing surfaces in the district itself
    for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j)
        textFile << "NA(NA):" << pDistrict->getSurface(j)->getId() << "(" << pDistrict->getSurface(j)->getKey() << "):" << unit << "\t";

    // loop on the obstructing surfaces in the district itself
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it)
        textFile << "NA(NA):" << (*it)->getId() << "(" << (*it)->getKey() << "):" << unit << "\t";

    textFile << endl;
    textFile.close();

}

void XmlScene::writeSWResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps (not giving the pre-time steps)
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        textFile << pClimate->getIgh(i) << "\t";
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                // loop on the walls, added the id of the surface
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getShortWaveIrradiance(i) << "\t";
                }
                // loop on the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getShortWaveIrradiance(i) << "\t";
                }
                // loop on the obstructing surfaces
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getShortWaveIrradiance(i) << "\t";
                }
            } // end the loop on zones
        }

        // loop on the tree surfaces in the district itself
        for (size_t j=0; j<pDistrict->getnTrees(); ++j)
            for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
                textFile << pDistrict->getTree(j)->getSurface(k)->getShortWaveIrradiance(i) << "\t";

        // loop on the obstructing surfaces in the district itself
        for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j)
            textFile << pDistrict->getSurface(j)->getShortWaveIrradiance(i) << "\t";

        // loop on the ground surfaces in the district itself
    for (forward_list<Ground*>::iterator it=getDistrict()->getGrounds()->begin();it!=getDistrict()->getGrounds()->end();++it)
            textFile << (*it)->getShortWaveIrradiance(i) << "\t";

        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeSWResultsBinary(string fileOut) {

    // the value as a float 32 bits
    float value;

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps (not giving the pre-time steps)
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        value = pClimate->getIgh(i);
        textFile.write(reinterpret_cast<char*>(&value),sizeof(value));
        // write the rest only if Igh > 0
        if (value > 0.f) {
            // loop on the number of buildings
            for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
                for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                    // loop on the walls, added the id of the surface
                    for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                        value = pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getShortWaveIrradiance(i);
                        textFile.write(reinterpret_cast<char*>(&value),sizeof(value));
                    }
                    // loop on the roofs
                    for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                        value = pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getShortWaveIrradiance(i);
                        textFile.write(reinterpret_cast<char*>(&value),sizeof(value));
                    }
                    // loop on the obstructing surfaces
                    for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                        value = pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getShortWaveIrradiance(i);
                        textFile.write(reinterpret_cast<char*>(&value),sizeof(value));
                    }
                } // end the loop on zones
            }
            // loop on the tree surfaces in the district itself
            for (size_t j=0; j<pDistrict->getnTrees(); ++j)
                for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k) {
                    value = pDistrict->getTree(j)->getSurface(k)->getShortWaveIrradiance(i);
                    textFile.write(reinterpret_cast<char*>(&value),sizeof(value));
                }
            // loop on the obstructing surfaces in the district itself
            for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j) {
                value = pDistrict->getSurface(j)->getShortWaveIrradiance(i);
                textFile.write(reinterpret_cast<char*>(&value),sizeof(value));
            }
            // loop on the ground surfaces in the district itself
            for (forward_list<Ground*>::iterator it=getDistrict()->getGrounds()->begin();it!=getDistrict()->getGrounds()->end();++it) {
                value = (*it)->getShortWaveIrradiance(i);
                textFile.write(reinterpret_cast<char*>(&value),sizeof(value));
            }
        }
    }
    textFile.close();

}

void XmlScene::writeSWvHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#" << "Igh_vis" << "\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << ":illuminance(lux)\t";
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getId() << ":illuminance(lux)\t";
            }
            // loop on the obstructing surfaces
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getId() << ":illuminance(lux)\t";
            }
        } // end the loop on zones
    }
    textFile << endl;
    textFile.close();

}

void XmlScene::writeSWvResultsText(string fileOut) {

    // test if something is to be written at all
    if (pClimate->isIgh_vis_empty()) return;

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps (not giving the pre-time steps)
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        textFile << pClimate->getIgh_vis(i) << "\t";
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                // loop on the walls, added the id of the surface
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getIlluminance(i) << "\t";
                }
                // loop on the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getIlluminance(i) << "\t";
                }
                // loop on the obstructing surfaces
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getIlluminance(i) << "\t";
                }
            } // end the loop on zones
        }
        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeDLHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#" << "Igh_vis" << "\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << ":internalIlluminanceNearWindow(lux)\t";
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << ":internalIlluminanceFarWindow(lux)\t";
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << ":internallyReflectedIlluminance(lux)\t";
            }
        } // end the loop on zones
    }
    textFile << endl;
    textFile.close();

}

void XmlScene::writeDLResultsText(string fileOut) {

    // test if something is to be written at all
    if (pClimate->isIgh_vis_empty()) return;

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        textFile << pClimate->getIgh_vis(i) << "\t";
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                // loop on the walls, added the id of the surface
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getInternalIlluminance0(i) << "\t";
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getInternalIlluminance1(i) << "\t";
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getInternallyReflectedIlluminance(i) << "\t";
                }
            } // end the loop on zones
        }
        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeLWHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << "\t";
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getId() << "\t";
            }
            // loop on the obstructing surfaces
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getId() << "\t";
            }
        } // end the loop on zones
    }

    // loop on the tree surfaces in the district itself
    for (size_t j=0; j<pDistrict->getnTrees(); ++j)
        for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
            textFile << pDistrict->getTree(j)->getId() << "(" << pDistrict->getTree(j)->getKey() << "):" << pDistrict->getTree(j)->getSurface(k)->getId() << "\t";

    // loop on the obstructing surfaces in the district itself
    for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j)
        textFile << "NA(NA):" << pDistrict->getSurface(j)->getId() << "\t";

    // loop on the ground surfaces in the district itself
    for (forward_list<Ground*>::iterator it=getDistrict()->getGrounds()->begin();it!=getDistrict()->getGrounds()->end();++it)
        textFile << "NA(NA):" << (*it)->getId() << "\t";

    textFile << endl;
    textFile.close();

}

void XmlScene::writeLWResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                // loop on the walls, added the id of the surface
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getLongWaveIrradiance(i) << "\t";
                }
                // loop on the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getLongWaveIrradiance(i) << "\t";
                }
                // loop on the obstructing surfaces
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getLongWaveIrradiance(i) << "\t";
                }
            } // end the loop on zones
        }

        // loop on the tree surfaces in the district itself
        for (size_t j=0; j<pDistrict->getnTrees(); ++j)
            for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
                textFile << pDistrict->getTree(j)->getSurface(k)->getLongWaveIrradiance(i) << "\t";

        // loop on the obstructing surfaces in the district itself
        for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j)
            textFile << pDistrict->getSurface(j)->getLongWaveIrradiance(i) << "\t";

        // loop on the obstructing surfaces in the district itself
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it)
            textFile << (*it)->getLongWaveIrradiance(i) << "\t";

        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeTHHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of buildings, to create the header
    textFile << "#timeStep\t";

    //beginning of contents added by Dapeng // Cognet: Adapted to new code.
    //loop on the number of DEC, to create the header
    for(unsigned int j=0; j<pDistrict->getnDECs(); ++j) {
        unsigned int decId = pDistrict->getDEC(j)->getId();
        textFile << "DEC" << decId <<":TotalThermalLoss(W)" <<"\t";
        pDistrict->getDEC(j)->getNetwork()->writeTHHeaderText(textFile);
    }
    //end of contents added by Dapeng

    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":Ta(celsius)" << "\t"
                     << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":Heating(Wh)" << "\t"
                     << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":Cooling(Wh)" << "\t";
            if (pDistrict->getBuilding(j)->getHVACpresence()) {
               textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACHeat" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACHeatAvailable" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACCool" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACCoolAvailable" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACHumidification" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACHumidificationAvailable" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACEvaporation" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACEvaporationAvailable" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACReheat" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACReheatAvailable" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACMassFlowRate" << "\t"
                        << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":HVACMassFlowRateAvailable" << "\t";
            }
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":Qi(Wh)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":Qs(Wh)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":VdotVent(m³/h)" << "\t";
        }
        textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):HeatStockTemperature(celsius)" << "\t";
        if (pDistrict->getBuilding(j)->getDHWHeatStock() != NULL)
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):DHWStockTemperature(celsius)" << "\t";
        textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):ColdStockTemperature(celsius)" << "\t"
                 << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):MachinePower(W)" << "\t"
                 << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):FuelConsumption(MJ)" << "\t"
                 << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):ElectricConsumption(kWh)" << "\t"
                 << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):SolarPVProduction(kWh)" << "\t"
                 << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):SolarThermalProduction(Wh)" << "\t";
    }
    textFile << endl;
    textFile.close();

}

void XmlScene::writeTHResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // if in pre-simulation phase, output all the values
    float preTimeStepsSimulated = this->preTimeStepsSimulated;
    float timeStepsSimulated = this->timeStepsSimulated;
    if (timeStepsSimulated == 0) { timeStepsSimulated = preTimeStepsSimulated; preTimeStepsSimulated = 0; }

    // loop on the number of time steps
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        textFile << simulationIndex+i+1-((preTimeStepsSimulated>0)?preTimeStepsSimulated:timeStepsSimulated) << "\t";

        // loop on the number of DEC, add by Dapeng // Cognet: Adapted to new code.
        for(unsigned int j=0; j<pDistrict->getnDECs(); ++j) {
            DistrictEnergyCenter* dec = pDistrict->getDEC(j);
            textFile << fixed << setprecision(0) << dec->getTotalThermalLoss(i)<<"\t";
            dec->getNetwork()->writeTHResultsText(textFile, i);
        }

        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            // loop on the zones
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                textFile << fixed << setprecision(1) << pDistrict->getBuilding(j)->getZone(k)->getTa(i) << "\t";
                textFile << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHeating(i-preTimeStepsSimulated+simulationIndex) << "\t";
                textFile << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getCooling(i-preTimeStepsSimulated+simulationIndex) << "\t";
                // if the HVAC system is present, save those values
                if (pDistrict->getBuilding(j)->getHVACpresence()) {
                   textFile << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACHeat(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACHeatAvailable(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACCool(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACCoolAvailable(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACHumidification(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACHumidificationAvailable(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACEvaporation(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACEvaporationAvailable(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACReheat(i) << "\t"
                            << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getHVACReheatAvailable(i) << "\t"
                            << fixed << setprecision(1) << pDistrict->getBuilding(j)->getZone(k)->getHVACMassFlowRate(i) << "\t"
                            << fixed << setprecision(1) << pDistrict->getBuilding(j)->getZone(k)->getHVACMassFlowRateAvailable(i) << "\t";
                }
                textFile << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getQi(i) << "\t";
                textFile << fixed << setprecision(0) << pDistrict->getBuilding(j)->getZone(k)->getQs(i-preTimeStepsSimulated) << "\t";
                textFile << fixed << setprecision(1) << 3600.f*pDistrict->getBuilding(j)->getZone(k)->getVdotVent(i) << "\t";
            }
            textFile << fixed << setprecision(1) << pDistrict->getBuilding(j)->getHeatStockTemperature(i) << "\t";
            if (pDistrict->getBuilding(j)->getDHWHeatStock() != NULL)
                textFile << fixed << setprecision(1) << pDistrict->getBuilding(j)->getDHWStockT(i) << "\t";
            textFile << fixed << setprecision(1) << pDistrict->getBuilding(j)->getColdStockTemperature(i) << "\t"
                     << fixed << setprecision(0) << pDistrict->getBuilding(j)->getMachinePower(i+simulationIndex) << "\t"
                     << fixed << setprecision(3) << pDistrict->getBuilding(j)->getFuelConsumption(i-preTimeStepsSimulated+simulationIndex)/1.e6 << "\t"
                     << fixed << setprecision(3) << pDistrict->getBuilding(j)->getElectricConsumption(i-preTimeStepsSimulated+simulationIndex)/3.6e6 << "\t"
                     << fixed << setprecision(0) << pDistrict->getBuilding(j)->getSolarPVProduction(i-preTimeStepsSimulated+simulationIndex)/3.6e6 << "\t"
                     << fixed << setprecision(0) << pDistrict->getBuilding(j)->getSolarThermalProduction(i-preTimeStepsSimulated+simulationIndex)/3.6e3 << "\t";

        }
        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeTHExplicitHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of buildings, to create the header
    textFile << "#timeStep\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":Ta(celsius)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":windowState" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":VdotVent(m3/s)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":lowerShadingState" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":upperShadingState" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(k)->getId() << ":lumint(lux)" << "\t";
        }
    }
    textFile << endl;
    textFile.close();

}

void XmlScene::writeTHExplicitResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps
    for (unsigned int i=preTimeSteps2Simulated; i<preTimeSteps2Simulated+timeSteps2Simulated; ++i) {
        textFile << static_cast<float>(simulationIndex+i-preTimeSteps2Simulated) / (Model::dt/Model::dt2) << "\t";
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            // loop on the zones
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                textFile << pDistrict->getBuilding(j)->getZone(k)->getTaExpl(i) << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(k)->getWindowState(i) << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(k)->getVdotVent(i) << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(k)->getTotalInternalIlluminance0(i) << "\t";
            }
        }
        textFile << endl;
    }

    textFile.close();

}

void XmlScene::writeTSHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // output of the Ke of that zone
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getId() << ":Ke(W/(m2K))\t";
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << ":Tos(°C)\t";
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getId() << ":Tos(°C)\t";
            }
        } // end the loop on zones
    }

    // loop on the tree surfaces in the district itself
    for (size_t j=0; j<pDistrict->getnTrees(); ++j)
        for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
            textFile << pDistrict->getTree(j)->getId() << "(" << pDistrict->getTree(j)->getKey() << "):" << pDistrict->getTree(j)->getSurface(k)->getId() << "\t";

    // loop on the ground surfaces
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        textFile << "NA(NA):" << (*it)->getId() << ":Tos(°C)\t";
    }

    textFile << endl;
    textFile.close();

}

void XmlScene::writeTSResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                // output of the Ke of that zone
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getKe(i)/pDistrict->getBuilding(j)->getZone(zone)->getSwa() << "\t";
                // loop on the walls, added the id of the surface
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) { // TODO: store the temperatures
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getTemperature(i) << "\t";
                }
                // loop on the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getTemperature(i) << "\t";
                }
            } // end the loop on zones
        }

        // loop on the tree surfaces in the district itself
        for (size_t j=0; j<pDistrict->getnTrees(); ++j)
            for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k)
                textFile << pDistrict->getTree(j)->getSurface(k)->getTemperature(i) << "\t";

        for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
            // loop on the ground surfaces
            textFile << (*it)->getTemperature(i) << "\t";
        }

        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeHCHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << ":hc(W/(m²K))\t";
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getId() << ":hc(W/(m²K))\t";
            }
        } // end the loop on zones
    }
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        // loop on the ground surfaces
        textFile << "NA(NA):" << (*it)->getId() << ":hc(W/(m²K))\t";
    }
    textFile << endl;
    textFile.close();

}

void XmlScene::writeHCResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                // loop on the walls, added the id of the surface
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->get_hc(i) << "\t";
                }
                // loop on the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                    textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->get_hc(i) << "\t";
                }
            } // end the loop on zones
        }
        for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
            // loop on the ground surfaces
            textFile << (*it)->get_hc(i) << "\t";
        }
        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeETHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                if (pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->hasET())
                    textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getId() << ":waterEvapotranspiration(mm/m²)\t";
            }
        } // end the loop on zones
    }
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        // loop on the ground surfaces
        if ((*it)->hasET())
            textFile << "NA(NA):" << (*it)->getId() << ":waterEvapotranspiration(mm/(m²h))\t";
    }
    textFile << endl;
    textFile.close();

}

void XmlScene::writeETResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop on the number of time steps
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
                // loop on the roofs
                for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                    if (pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->hasET())
                        textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getWaterEvapotranspiration(i) << "\t";
                }
            }
        }
        // loop on the ground surfaces
        for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
            if ((*it)->hasET())
                textFile << (*it)->getWaterEvapotranspiration(i) << "\t";
        }
        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeCMHeaderText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // loop to output the surface ids
    textFile << "#";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            if (!pDistrict->getBuilding(j)->isMRT()) continue; // forgets about the buildings that are not MRT (goes at the end of the for statement)
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << ")" << ":MRT(°C)\t"
                     << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << ")" << ":COMFA*(W/m²)\t"
                     << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << ")" << ":ITS(W)\t";
    }

    textFile << endl;
    textFile.close();

}

void XmlScene::writeCMResultsText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);
    if (!textFile.is_open()) throw(string("Cannot open file: ")+fileOut);

    // variables used
    float MRT, COMFA, ITS;

    // loop on the number of time steps
    for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            if (!pDistrict->getBuilding(j)->isMRT()) continue; // forgets about the buildings that are not MRT (goes at the end of the for statement)

            // calls the model to compute the comfort indices
            Model::computeCMIndices(pClimate,pDistrict->getBuilding(j),i,preTimeStepsSimulated,MRT,COMFA,ITS);

            textFile << MRT << "\t";
            textFile << COMFA << "\t";
            textFile << ITS << "\t";

        }
        textFile << endl;
    }
    textFile.close();

}

void XmlScene::writeClimaticDataText(string fileOut) {

    // clear the file
    remove(fileOut.c_str());

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);

    // write the header
    textFile << "#day\thour\tTout\tTd\tTg\tTsky\tPa\tg_vs_ab" << endl;

    // debut de la simulation sur la periode consideree
    for (unsigned int day = beginDay; day<=endDay; ++day) {
        for (unsigned int hour = 1; hour <= 24; ++hour) {
            textFile << day << "\t" << hour << "\t" << pClimate->getToutCelsius(day,hour) << "\t" << pClimate->getTd(day,hour) << "\t"
                     << pClimate->getTgroundCelsius(day,hour) << "\t" << pClimate->getTskyCelsius(day,hour) << "\t" << pClimate->getPatm(day,hour)/1000.
                     << "\t" << toString(((pClimate->getIdh(day,hour)>0.f)||(pClimate->getRelativeHumidity(day,hour)>0.8f))?0.3f:0.01f) << endl;
        }
    }

    textFile.close();

}

void XmlScene::writeMonthlyResultsText(string fileOut) {

    // clear the file
    remove(fileOut.c_str());

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);

    int daysPerMonth[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    int cumDays[] = {0,31,59,90,120,151,181,212,243,273,304,334,365};

    // loop on the number of buildings, to create the header
    textFile << "#month\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << "heatingNeeds(Wh)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << "coolingNeeds(Wh)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << "DHW(J)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << "fuelConsumption(J)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << "electricConsumption(J)" << "\t";
            textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << "electricProduction(J)" << "\t";
    }
    textFile << endl;

    for (unsigned int monthIndex=0; monthIndex<12; ++monthIndex) {

        textFile << monthIndex+1 << "\t";
        // loop on the number of buildings
        for (unsigned int j=0; j<pDistrict->getnBuildings(); j++) {
            float heating = 0.f;
            float cooling = 0.f;
            float dhw = 0.f;
            float fuel = 0.f;
            float electric = 0.f;
            float pv = 0.f;
            dhw = 1.2f*1000.f*50.e-3f*daysPerMonth[monthIndex]*pDistrict->getBuilding(j)->getOccupantsCount()*(55.f-10.f);
            // loop on the number of time steps
            for (unsigned int i=preTimeStepsSimulated+cumDays[monthIndex]*24; i<preTimeStepsSimulated+cumDays[monthIndex+1]*24; ++i) {
                heating += pDistrict->getBuilding(j)->getHeating(i-preTimeStepsSimulated);
                cooling += pDistrict->getBuilding(j)->getCooling(i-preTimeStepsSimulated);
                fuel += pDistrict->getBuilding(j)->getFuelConsumption(i-preTimeStepsSimulated);
                electric += pDistrict->getBuilding(j)->getElectricConsumption(i-preTimeStepsSimulated);
                pv += pDistrict->getBuilding(j)->getSolarPVProduction(i-preTimeStepsSimulated);
            }
            textFile << heating << "\t" << cooling << "\t" << dhw << "\t" << fuel << "\t" << electric << "\t";
        }
        textFile << endl;
    }

    textFile.close();

}

void XmlScene::writeYearlyResultsPerBuildingText(string fileOut) {

    // clear the file
    remove(fileOut.c_str());

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);

    // write the header
    textFile << "#buildingId(key)\theatingNeeds(Wh)\tcoolingNeeds(Wh)" << endl;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        double heating = 0.;
        double cooling = 0.;
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            heating += pDistrict->getBuilding(j)->getZone(k)->getTotalHeatingSatisfied();
            cooling += pDistrict->getBuilding(j)->getZone(k)->getTotalCoolingSatisfied();
        }
        textFile << pDistrict->getBuilding(j)->getId()
                 << "(" << pDistrict->getBuilding(j)->getKey() << "):"
                 << "\t" << heating << "\t" << cooling << endl;
    }
    textFile.close();

}

////added by Dapeng // Cognet: Deleted it (never used?).
//void XmlScene::writeYearlyResultsPerDECText(string fileOut) {
//    // clear the fiel
//    remove(fileOut.c_str());

//    //open the binary file
//    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::app);

//    //write the header
//    textFile<<"#DECID\t fuelConsumption\t electricConsumption"<<endl;

//    //loop on the number of DECs
//    for (unsigned int j=0; j<pDistrict->getnDECs(); ++j) {
//        textFile<<pDistrict->getDEC(j)->getId()<<"\t"
//                <<pDistrict->getDEC(j)->getTotalFuelConsumption()<<"\t"
//                <<pDistrict->getDEC(j)->getTotalElectricConsumption()<<endl;
//    }
//    textFile.close();

//}

map<unsigned int,vector<double> > XmlScene::getHeatingHourlyResultsPerBuilding() {

    map<unsigned int,vector<double> > results;
    vector<double> heating;
    double heatingHour_i = 0.;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        heating.clear();
        // loop on the simulation time steps
        for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+simulationIndex+timeStepsSimulated; ++i) { // DP: the total length of the vector is simulationIndex + timeStepsSimulated
            heatingHour_i = 0.;
            // loop on the zones
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                heatingHour_i += pDistrict->getBuilding(j)->getZone(k)->getHeating(i-preTimeStepsSimulated); // DP: Heating does not keep any preTimeStepsSimulated values
            }
            heating.push_back(heatingHour_i);
        }
        results.insert(pair<unsigned int,vector<double> >(pDistrict->getBuilding(j)->getId(), heating));
    }

    return results;

}

map<unsigned int,vector<double> > XmlScene::getHeatingMonthlyResultsPerBuilding() {

    map<unsigned int,vector<double> > results;

    // days in months
    int cumDays[] = {0,31,59,90,120,151,181,212,243,273,304,334,365};

    for (unsigned int j=0; j<pDistrict->getnBuildings(); j++) {

        vector<double> heating;
        // loop on the number of buildings
        for (unsigned int monthIndex=0; monthIndex<12; ++monthIndex) {
            double heatingPerMonth = 0.;
            // loop on the zones per building
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                // loop on the simulation time steps
                for (unsigned int i=preTimeStepsSimulated+cumDays[monthIndex]*24; i<preTimeStepsSimulated+cumDays[monthIndex+1]*24; ++i) {
                    heatingPerMonth += pDistrict->getBuilding(j)->getZone(k)->getHeating(i-preTimeStepsSimulated); // DP: Heating does not keep any preTimeStepsSimulated values
                }
            }
            heating.push_back(heatingPerMonth);
        }
        // puts back the result in the set
        results.insert(pair<unsigned int,vector<double> >(pDistrict->getBuilding(j)->getId(), heating));
    }

    return results;

}

map<unsigned int,double> XmlScene::getHeatingYearlyResultsPerBuilding() {

    map<unsigned int,double> results;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        double heating = 0.;
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            heating += pDistrict->getBuilding(j)->getZone(k)->getTotalHeating();
        }
        results.insert(pair<unsigned int,double>(pDistrict->getBuilding(j)->getId(), heating));
    }

    return results;

}

map<unsigned int,vector<double> > XmlScene::getCoolingHourlyResultsPerBuilding() {

    map<unsigned int,vector<double> > results;
    vector<double> cooling;
    double coolingHour_i = 0.;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        cooling.clear();
        // loop on the simulation time steps
        for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+simulationIndex+timeStepsSimulated; ++i) { // DP: the total length of the vector is simulationIndex + timeStepsSimulated
            coolingHour_i = 0.;
            // loop on the zones
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                coolingHour_i += pDistrict->getBuilding(j)->getZone(k)->getCooling(i-preTimeStepsSimulated); // DP: Heating does not keep any preTimeStepsSimulated values
            }
            cooling.push_back(coolingHour_i);
        }
        results.insert(pair<unsigned int,vector<double> >(pDistrict->getBuilding(j)->getId(), cooling));
    }

    return results;

}

map<unsigned int,vector<double> > XmlScene::getCoolingMonthlyResultsPerBuilding() {

    map<unsigned int,vector<double> > results;

    // days in months
    int cumDays[] = {0,31,59,90,120,151,181,212,243,273,304,334,365};

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); j++) {

        vector<double> cooling;
        // loop on months
        for (unsigned int monthIndex=0; monthIndex<12; ++monthIndex) {
            double coolingPerMonth = 0.;
            // loop on the zones per building
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
                // loop on the number of time steps
                for (unsigned int i=preTimeStepsSimulated+cumDays[monthIndex]*24; i<preTimeStepsSimulated+cumDays[monthIndex+1]*24; ++i) {
                    coolingPerMonth += pDistrict->getBuilding(j)->getZone(k)->getCooling(i-preTimeStepsSimulated); // DP: Heating does not keep any preTimeStepsSimulated values
                }
            }
            cooling.push_back(coolingPerMonth);
        }
        // puts back the result in the set
        results.insert(pair<unsigned int,vector<double> >(pDistrict->getBuilding(j)->getId(), cooling));
    }

    return results;

}

map<unsigned int,double> XmlScene::getCoolingYearlyResultsPerBuilding() {

    map<unsigned int,double> results;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        double cooling = 0.;
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            cooling += pDistrict->getBuilding(j)->getZone(k)->getTotalCooling();
        }
        results.insert(pair<unsigned int,double>(pDistrict->getBuilding(j)->getId(), cooling));
    }

    return results;

}

map<unsigned int,vector<double> > XmlScene::getMachinePowerHourlyResultsPerBuilding() {

    map<unsigned int,vector<double> > results;
    vector<double> qs;
    //double qsHour_i = 0.;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        qs.clear();
        // loop on the time steps
        for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+simulationIndex+timeStepsSimulated; ++i) { // DP: the total length of the vector is preTimeStepsSimulated + simulationIndex + timeStepsSimulated
            qs.push_back(pDistrict->getBuilding(j)->getMachinePower(i)); // DP: MachinePower contains preTimeStepsSimulated=1 pre-simulated value at the begining
        }
        results.insert(pair<unsigned int,vector<double> >(pDistrict->getBuilding(j)->getId(), qs));
    }

    return results;

}

map<unsigned int,vector<double> > XmlScene::getMachinePowerMonthlyResultsPerBuilding() {

    map<unsigned int,vector<double> > results;

    // days in months
    int cumDays[] = {0,31,59,90,120,151,181,212,243,273,304,334,365};

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); j++) {

        vector<double> qs;
        // loop on months
        for (unsigned int monthIndex=0; monthIndex<12; ++monthIndex) {
            double qsPerMonth = 0.;
            // loop on the number of time steps
            for (unsigned int i=preTimeStepsSimulated+cumDays[monthIndex]*24; i<preTimeStepsSimulated+cumDays[monthIndex+1]*24; ++i) { // DP: the total length of the vector is preTimeStepsSimulated + simulationIndex + timeStepsSimulated
                qsPerMonth += pDistrict->getBuilding(j)->getMachinePower(i); // DP: MachinePower contains preTimeStepsSimulated=1 pre-simulated value at the begining
            }
            qs.push_back(qsPerMonth);
        }
        // puts back the result in the set
        results.insert(pair<unsigned int,vector<double> >(pDistrict->getBuilding(j)->getId(), qs));
    }
    return results;

}

map<unsigned int,double> XmlScene::getMachinePowerYearlyResultsPerBuilding() {

    map<unsigned int,double> results;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); j++) {
        double qs=0;

        // loop on the number of time steps
        for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+simulationIndex+timeStepsSimulated; ++i) { // DP: the total length of the vector is preTimeStepsSimulated + simulationIndex + timeStepsSimulated
            qs += pDistrict->getBuilding(j)->getMachinePower(i); // DP: MachinePower contains preTimeStepsSimulateds=1 pre-simulated value at the begining
        }
        // puts back the result in the set
        results.insert(pair<unsigned int,double>(pDistrict->getBuilding(j)->getId(), qs));
    }
    return results;

}

void XmlScene::writeYearlyResultsText(string fileOut) {

    // open the binary file truncating the content of the previous file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);

    float heating = 0.f;
    float cooling = 0.f;
    float dhw = 0.f;
    float fuel = 0.f;
    float electric = 0.f;
    float pv = 0.f;
    // eco indicators
    float nre = 0.f;
    float gwp = 0.f;
    float ubp = 0.f;
    // MRT indicator
    float mrt = 0.f, MRTcount=0.f;
    float MRT, COMFA, ITS, COMFAcount=0.f, ITScount=0.f;
    // KPI
    float totalThermalLoss = 0.f;
    float totalHeatingUnsatisfied = 0.f;
    float totalCoolingUnsatisfied = 0.f;
    float electricPump = 0.f;

    // loop on the number of buildings
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        dhw += 1.2f*1000.f*50.e-3f*365.f*pDistrict->getBuilding(j)->getOccupantsCount()*(55.f-10.f);
        // eco indicators
        nre += pDistrict->getBuilding(j)->getNRE();
        gwp += pDistrict->getBuilding(j)->getGWP();
        ubp += pDistrict->getBuilding(j)->getUBP();
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            heating += pDistrict->getBuilding(j)->getZone(k)->getTotalHeatingSatisfied();
            cooling += pDistrict->getBuilding(j)->getZone(k)->getTotalCoolingSatisfied();
        }
        fuel += pDistrict->getBuilding(j)->getTotalFuelConsumption();
        electric += pDistrict->getBuilding(j)->getTotalElectricConsumption();
        pv += pDistrict->getBuilding(j)->getTotalSolarPVProduction();
        // MRT indicator
        for (unsigned int i=preTimeStepsSimulated; i<preTimeStepsSimulated+timeStepsSimulated; ++i) {
            if (!pDistrict->getBuilding(j)->isMRT()) continue; // forgets about the buildings that are not MRT (goes at the end of the for statement)

            // compute the comfort indices
            Model::computeCMIndices(pClimate,pDistrict->getBuilding(j),i,preTimeStepsSimulated,MRT,COMFA,ITS);
            // keep only the values between 8h and 18h
            if  ( ((i-preTimeStepsSimulated)%24 >= 7) && ((i-preTimeStepsSimulated)%24 <= 17) ) {
                // computes the average mrt for the defined hours
                mrt += MRT;
                MRTcount += 1.f;
                // computes the number of hours of comfort for ITS
                if (ITS > -160.f && ITS < 160.f)
                    ITScount += 1.f;
                // computes the number of hours of comfort for COMFA
                if (COMFA > -50.f && COMFA < 50.f)
                    COMFAcount += 1.f;

            }
        }
        // KPI
        totalHeatingUnsatisfied += pDistrict->getBuilding(j)->getHeatingDemandUnsatisfied();
        totalCoolingUnsatisfied += pDistrict->getBuilding(j)->getCoolingDemandUnsatisfied();
    }

    //Added by Max
    for(size_t i(0); i < pDistrict->getnDECs(); ++i){
        totalThermalLoss += pDistrict->getDEC(i)->getYearlyTotalThermalLoss();
        for(size_t j(0); j < pDistrict->getDEC(i)->getnThermalStations(); ++j){
           electricPump += pDistrict->getDEC(i)->getThermalStation(j)->getTotalPumpPower();
        }
    }

    textFile << "#Total yearly energy demand & supply\n";
    textFile << "Heating (Wh):\t" << heating << "\n";
    textFile << "Cooling (Wh):\t" << cooling << "\n";
    textFile << "DHW (J):\t" << dhw << "\n";
    textFile << "Fuel (J):\t" << fuel << "\n";
    textFile << "ElectricConsumption (J):\t" << electric << "\n";
    textFile << "ElectricPVProduction (J):\t" << pv << endl;
    textFile << "HeatingUnsatisfied (Wh):\t" << totalHeatingUnsatisfied << endl; //Added by Max
    textFile << "CoolingUnsatisfied (Wh):\t" << totalCoolingUnsatisfied << endl; //Added by Max
    textFile << "NRE (MJ):\t" << nre << "\n";
    textFile << "GWP (kgCO2):\t" << gwp << "\n";
    textFile << "UBP (pts):\t" << ubp << endl;
    textFile << "MRT (°C):\t" << mrt/MRTcount << endl;
    textFile << "ITS (h):\t" << ITScount << endl;
    textFile << "COMFA (h):\t" << COMFAcount << endl;
    textFile << "DHNTotalThermalLosses (Wh):\t" << totalThermalLoss << endl; //Added by Max
    textFile << "DHNElectricPump (J):\t" << electricPump << endl; //Added by Max

// Cognet: Deleted it, to put back, need to change how Fuel and ElectricConsumption are saved.
//    //beginning of contents added by Dapeng
//    double fuel_DEC=0, electric_DEC=0;
//    for(unsigned int i=0; i<pDistrict->getnDECs(); ++i) {
//        fuel_DEC += pDistrict->getDEC(i)->getTotalFuelConsumption();
//        electric_DEC += pDistrict->getDEC(i)->getTotalElectricConsumption();
//    }
//    if(pDistrict->getnDECs() > 0) {
//        textFile<<"District Fuel (J):\t" << fuel_DEC <<"\n";
//        textFile<<"District Electric (J):\t" << electric_DEC <<"\n";
//    }
//    //end of contents added by Dapeng

//    logStream << "Total yearly energy demand & supply" << endl << flush;
//    logStream << "heating (Wh): " << heating << "\ncooling (Wh): " << cooling << "\nfuel (J): " << fuel
//         << "\ndhw (J): " << dhw << "\nelectric (J): " << electric << endl << flush;

    textFile.close();

}

void XmlScene::writeAreaText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);

    // loop on the number of buildings, to create the header
    textFile << "type\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << ")\t";
    }
    textFile << endl;

    // loop on the number of buildings, to show the result for the FLOOR area
    textFile << "floorArea(m2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getFloorArea() << "\t";
    }
    textFile << endl;
    // NRE
    textFile << "floorNRE(MJ)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getFloorNRE() << "\t";
    }
    textFile << endl;
    // GWP
    textFile << "floorGWP(kgCO2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getFloorGWP() << "\t";
    }
    textFile << endl;
    // UBP
    textFile << "floorUBP(points)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getFloorUBP() << "\t";
    }
    textFile << endl;

    // loop on the number of buildings, to show the result for the ROOF area
    textFile << "roofArea(m2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getRoofArea() << "\t";
    }
    textFile << endl;
    // PV
    textFile << "roofPVArea(m2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getRoofPVArea() << "\t";
    }
    textFile << endl;
    // NRE
    textFile << "roofNRE(MJ)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getRoofNRE() << "\t";
    }
    textFile << endl;
    // GWP
    textFile << "roofGWP(kgCO2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getRoofGWP() << "\t";
    }
    textFile << endl;
    // UBP
    textFile << "roofUBP(points)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getRoofUBP() << "\t";
    }
    textFile << endl;

    // loop on the number of buildings, to show the result for the wall area
    textFile << "wallArea(m2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getWallArea() << "\t";
    }
    textFile << endl;
    // PV
    textFile << "wallPVArea(m2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getWallPVArea() << "\t";
    }
    textFile << endl;
    // NRE
    textFile << "wallNRE(MJ)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getWallNRE() << "\t";
    }
    textFile << endl;
    // GWP
    textFile << "wallGWP(kgCO2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getWallGWP() << "\t";
    }
    textFile << endl;
    // UBP
    textFile << "wallUBP(points)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getWallUBP() << "\t";
    }
    textFile << endl;

    // loop on the number of buildings, to show the result for the window areas
    textFile << "windowsArea(m2)\t";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getWindowsArea() << "\t";
    }
    textFile << endl;

    textFile.close();

}

void XmlScene::writeVFText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);

    // 0. writes the header of the file
    textFile << "#\tType\tArea(m²)\tAzimuth(°)\tAltitude(°)\tSVF(-)\tGVF(-)" << endl;

    // 1. write for each surface the information in the header
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        for (unsigned int zone=0; zone<pDistrict->getBuilding(j)->getnZones();++zone) {
            // loop on the walls, added the id of the surface
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnWalls(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getId() << "\t";
                textFile << "Wall" << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getArea() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getAzimuth() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getAltitude() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getSVF() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getWall(k)->getGVF() << endl;
            }
            // loop on the roofs
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnRoofs(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getId() << "\t";
                textFile << "Roof" << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getArea() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getAzimuth() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getAltitude() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getSVF() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getRoof(k)->getGVF() << endl;

            }
            // loop on the obstructing surfaces
            for (unsigned int k=0; k<pDistrict->getBuilding(j)->getZone(zone)->getnSurfaces(); ++k) {
                textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << "):" << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getId() << "\t";
                textFile << "Surface" << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getArea() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getAzimuth() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getAltitude() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getSVF() << "\t";
                textFile << pDistrict->getBuilding(j)->getZone(zone)->getSurface(k)->getGVF() << endl;
            }
        } // end the loop on zones
    }
    // loop on the tree surfaces in the district itself
    for (size_t j=0; j<pDistrict->getnTrees(); ++j) {
        for (size_t k=0;k<pDistrict->getTree(j)->getnSurfaces();++k) {
            textFile << pDistrict->getTree(j)->getId() << "(" << pDistrict->getTree(j)->getKey() << "):" << pDistrict->getTree(j)->getSurface(k)->getId() << "\t";
            textFile << "Tree" << "\t";
            textFile << pDistrict->getTree(j)->getSurface(k)->getArea() << "\t";
            textFile << pDistrict->getTree(j)->getSurface(k)->getAzimuth() << "\t";
            textFile << pDistrict->getTree(j)->getSurface(k)->getAltitude() << "\t";
            textFile << pDistrict->getTree(j)->getSurface(k)->getSVF() << "\t";
            textFile << pDistrict->getTree(j)->getSurface(k)->getGVF() << endl;
        }
    }

    // loop on the obstructing surfaces in the district itself
    for (unsigned int j=0; j<pDistrict->getnSurfaces(); ++j) {
        textFile << "NA(NA):" << pDistrict->getSurface(j)->getId() << "\t";
        textFile << "Surface" << "\t";
        textFile << pDistrict->getSurface(j)->getArea() << "\t";
        textFile << pDistrict->getSurface(j)->getAzimuth() << "\t";
        textFile << pDistrict->getSurface(j)->getAltitude() << "\t";
        textFile << pDistrict->getSurface(j)->getSVF() << "\t";
        textFile << pDistrict->getSurface(j)->getGVF() << endl;
    }

    // loop on the obstructing surfaces in the district itself
    for (forward_list<Ground*>::iterator it=pDistrict->getGrounds()->begin();it!=pDistrict->getGrounds()->end();++it) {
        textFile << "NA(NA):" << (*it)->getId() << "\t";
        textFile << "Ground" << "\t";
        textFile << (*it)->getArea() << "\t";
        textFile << (*it)->getAzimuth() << "\t";
        textFile << (*it)->getAltitude() << "\t";
        textFile << (*it)->getSVF() << "\t";
        textFile << (*it)->getGVF() << endl;
    }

    textFile << endl;
    textFile.close();

}

void XmlScene::writeInertiaText(string fileOut) {

    // open the binary file
    fstream textFile(fileOut.c_str(),ios::out|ios::binary|ios::trunc);

    textFile << "#building\twarmUpTime(days)\tzonesThermalInertia(J/K)\n";
    for (unsigned int j=0; j<pDistrict->getnBuildings(); ++j) {
        textFile << pDistrict->getBuilding(j)->getId() << "(" << pDistrict->getBuilding(j)->getKey() << ")\t" << Model::ThermalWarmUpTime(pDistrict->getBuilding(j)) << "\t";
        // loop on the zones
        for (unsigned int k=0; k<pDistrict->getBuilding(j)->getnZones(); ++k) {
            for (unsigned int l=0; l<pDistrict->getBuilding(j)->getZone(k)->getnNodes(); ++l)
                textFile << pDistrict->getBuilding(j)->getZone(k)->getC(l) << "\t";
        }
        textFile << endl;
    }

    textFile.close();

}

void XmlScene::exportCumulativeRadiance()
{
    // computes the sky and ground file
    computeCumulativeRadiance(beginDay,endDay);
    computeFarField();
    exportSkyAndGround((inputFile.substr(0,inputFile.size()-4) + "_cumSkyGrnd.rad"));

    // preparation of the cumulative suns file
    ofstream radOut((inputFile.substr(0,inputFile.size()-4) + "_cumSuns.rad").c_str(), ios::binary);
    if (!radOut.is_open()) throw string("Error creating the cumulative Suns file.");

    double sunLuminance = 0.;
    GENPoint sunPosition;
    for (unsigned int day = beginDay; day<=endDay; ++day) {
        for (unsigned int hour = 1; hour <= 24; ++hour) {

            // initialisation of the sun (for the VFC)
            pSun->SetDay(day);
            pSun->SetClockTime1(hour);
            if (pSun->SunUp() && sunVisibleFarField(pSun->GetPosition().Azimuth().degrees(), pSun->GetPosition().Altitude().degrees())) { // if sun is up and visible

                sunLuminance = pClimate->getIbn(day,hour) / pSun->getSolidAngle();

                radOut << "void light solar_" << day << "_" << hour << "\n"
                       << "0\n"
                       << "0\n"
                       << "3 " << sunLuminance << " " << sunLuminance << " " << sunLuminance << endl;

                sunPosition = pSun->GetPosition();

                radOut << "solar_" << day << "_" << hour << " source sun_" << day << "_" << hour << "\n"
                       << "0\n"
                       << "0\n"
                       << "4 " << sunPosition[0] << " "
                               << sunPosition[1] << " "
                               << sunPosition[2] << " "
                               << pSun->getAperture() << "\n"
                       << endl;
            }
        }
    }

}

void XmlScene::exportHourlyRadiance() {

    // N.B.: the far Field obstructions vector needs to have been initialized before

    double sunLuminance = 0.;
    GENPoint sunPosition;
    for (unsigned int day = beginDay; day<=endDay; ++day) {
        for (unsigned int hour = 1; hour <= 24; ++hour) {

            // computes the radiance distribution at that time
            computeRadiance(day,pClimate->getIdh(day,hour),pClimate->getIbn(day,hour),pDistrict->getGroundAlbedo());
            computeFarField();
            exportSkyAndGround(inputFile.substr(0,inputFile.size()-4) + "_" + to_string(day) + "_" + to_string(hour)  + ".rad");

            // preparation of the cumulative suns file
            ofstream radOut((inputFile.substr(0,inputFile.size()-4) + "_" + to_string(day) + "_" + to_string(hour) + ".rad").c_str(), ios::app | ios::binary);
            if (!radOut.is_open()) throw string("Error creating the sun file on " + to_string(day) + "_" + to_string(hour));
            radOut << endl;

            // initialisation of the sun (for the VFC)
            pSun->SetDay(day);
            pSun->SetClockTime1(hour);
            if (pSun->SunUp() && sunVisibleFarField(pSun->GetPosition().Azimuth().degrees(), pSun->GetPosition().Altitude().degrees())) { // if sun is up and visible

                sunLuminance = pClimate->getIbn(day,hour) / pSun->getSolidAngle();

                radOut << "\nvoid light solar_" << day << "_" << hour << "\n"
                       << "0\n"
                       << "0\n"
                       << "3 " << sunLuminance << " " << sunLuminance << " " << sunLuminance << endl;

                sunPosition = pSun->GetPosition();

                radOut << "\nsolar_" << day << "_" << hour << " source sun_" << day << "_" << hour << "\n"
                       << "0\n"
                       << "0\n"
                       << "4 " << sunPosition[0] << " "
                               << sunPosition[1] << " "
                               << sunPosition[2] << " "
                               << pSun->getAperture() << "\n"
                       << endl;
            }
            radOut.close();
        }
    }

}

void XmlScene::exportSkyRadFile(string skyRadFile, float *lv) {

    // create the rad file with reference to .cal file

    ofstream skyRadOut(skyRadFile.c_str(), ios::binary);
    if (!skyRadOut.is_open()) throw(string("Error creating sky file: " + skyRadFile));

    skyRadOut << "#Sky file generated by CitySim (jerome.kaempf@kaemco.ch)"
              << "\n" << endl;

    skyRadOut << "void brightfunc skyfunc" << endl;
    skyRadOut << "2 skybright " << skyRadFile.substr(0,skyRadFile.size()-4) << ".cal" << endl;
    skyRadOut << "0" << endl << "0" << endl;
    skyRadOut << "skyfunc glow sky_glow" << endl;
    skyRadOut << "0" << endl << "0" << endl << "4 1 1 1 0" << endl;
    skyRadOut << "sky_glow source sky" << endl;
    skyRadOut << "0" << endl << "0" << endl << "4 0 0 1 180" << endl;

    skyRadOut.close();

    // create the .cal file with the luminances values

    // output the radiances
    ofstream skyCalOut((skyRadFile.substr(0,skyRadFile.size()-4) + ".cal").c_str(), ios::binary);
    if (!skyCalOut.is_open()) throw(string("Error creating sky file: " + skyRadFile.substr(0,skyRadFile.size()-4) + ".cal"));

    skyCalOut << "skybright=";
    for (unsigned int j=0; j<tregenzaSky.getBands()-1; j++)
    {
        skyCalOut << "row" << j << "+";
    }
    skyCalOut << "row" << tregenzaSky.getBands()-1 << ";" << endl << endl;

    unsigned int counter = 0;
    for (unsigned int j=0; j<tregenzaSky.getBands()-1; j++)
    {
        // note first patch split into two parts - first part (> 0 deg) and last patch (<360)
        skyCalOut << "row" << j << "=if(and(alt-" << j*tregenzaSky.getDeltaAltitude()*180./M_PI << ", " << (j+1)*tregenzaSky.getDeltaAltitude()*180./M_PI << "-alt),";
        skyCalOut << "select(floor(0.5+az/" << tregenzaSky.getDeltaAzimuth(j)*180./M_PI << ")+1," << endl;

        for (unsigned int i=counter; i< counter + tregenzaSky.getPatchesPerBand(j); i++)
        {
            skyCalOut << "\t" << lv[i] << "," << endl;
        }
        // rewrite the first one.
        skyCalOut << "\t" << lv[counter] << "),0);" << endl << endl;
        counter += tregenzaSky.getPatchesPerBand(j);
    }

    // top patch.
    skyCalOut << "row" << tregenzaSky.getBands()-1 << "=if(alt-"<< 90.-(tregenzaSky.getDeltaAltitude()*180./M_PI / 2.)<< "," << lv[counter] << ",0);"<< endl << endl;

    skyCalOut << "alt=asin(Dz)*180/PI;" << endl << endl;
    skyCalOut << "az=if(azi,azi,azi+360);" << endl;
    skyCalOut << "azi=atan2(Dx,Dy)*180/PI;" << endl << endl;

    skyCalOut.close();

}

size_t XmlScene::memoryUsage() {                    // DP: modified to take into account vectors that are not emptied
    // add the obstructing surfaces' 4 vectors of float (SW)
    size_t bytes = pDistrict->getnSurfaces()*sizeof(float)*1*(timeStepsSimulated+preTimeStepsSimulated);
    // add the ground surfaces' 4 vectors of float (SW,DLout,LW,TS)
    bytes += pDistrict->getnGrounds()*sizeof(float)*4*(timeStepsSimulated+preTimeStepsSimulated);
    // add the tree surfaces results
    for (size_t treeIndex=0; treeIndex<pDistrict->getnTrees(); ++treeIndex)
        bytes += pDistrict->getTree(treeIndex)->getnSurfaces()*sizeof(float)*3*(timeStepsSimulated+preTimeStepsSimulated); // TS, SW and LW
    // add for the buildings NOT EMPTIED vectors of double (machineP, el+fuelConsump)
    bytes += pDistrict->getnBuildings()*sizeof(double)*3*(timeStepsSimulated+preTimeStepsSimulated+simulationIndex);
    // add for the buildings EMPTIED vectors of double (2x stockT)
    bytes += pDistrict->getnBuildings()*sizeof(double)*2*(timeStepsSimulated+preTimeStepsSimulated);
    for (size_t buildingIndex=0; buildingIndex<pDistrict->getnBuildings(); ++buildingIndex) {
        // for each zone sum of NOT EMPTIED vectors of double (heating, cooling)
        bytes += pDistrict->getBuilding(buildingIndex)->getnZones()*sizeof(double)*2*(timeStepsSimulated+preTimeStepsSimulated+simulationIndex);
        // for each zone sum of EMPTIED vectors of double (Qs)
        bytes += pDistrict->getBuilding(buildingIndex)->getnZones()*sizeof(double)*1*(timeStepsSimulated+preTimeStepsSimulated);
        // DP: all other vectors are emptied and thus of size (timeStepsSimulated+preTimeStepsSimulated)
        // if HVAC present, some more results
        if (pDistrict->getBuilding(buildingIndex)->getHVACpresence()) {
            bytes += pDistrict->getBuilding(buildingIndex)->getnZones()*sizeof(double)*12*(timeStepsSimulated+preTimeStepsSimulated);
        }
        // for each zone sum of vectors of float (VdotVent,LowerShadingState,UpperShadingState,Lumint)
        if (Model::thermalExplicit)
            bytes += pDistrict->getBuilding(buildingIndex)->getnZones()*sizeof(float)*5*(timeSteps2Simulated+preTimeSteps2Simulated);
        for (size_t zoneIndex=0; zoneIndex<pDistrict->getBuilding(buildingIndex)->getnZones(); ++zoneIndex) {
            // for each zone sum of vectors for each node (Ta - 1N,Ta Tw - 2N or Ta Tw Tw2 - 3N) - Normal (hourly) + Explicit (5 mins)
            bytes += sizeof(double)*pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getnNodes()*(timeStepsSimulated+timeSteps2Simulated+preTimeStepsSimulated+preTimeSteps2Simulated);
            for (size_t wallIndex=0; wallIndex<pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getnWalls(); ++wallIndex) {
                // for each wall sum of 7 vectors of float (SW,DLout,LW,3*DL,TS)
                bytes += sizeof(*(pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getWall(wallIndex)));
                bytes += sizeof(float)*7*(timeStepsSimulated+preTimeStepsSimulated);
            }
            for (size_t roofIndex=0; roofIndex<pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getnRoofs(); ++roofIndex) {
                // for each wall sum of 4 vectors of float (SW,DLout,LW,TS)
                bytes += sizeof(*(pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getRoof(roofIndex)));
                bytes += sizeof(float)*4*(timeStepsSimulated+preTimeStepsSimulated);
            }
            for (size_t surfaceIndex=0; surfaceIndex<pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getnSurfaces(); ++surfaceIndex) {
                // for each wall sum of 3 vectors of float (SW,DLout,LW) and no TS
                bytes += sizeof(*(pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getSurface(surfaceIndex)));
                bytes += sizeof(float)*3*(timeStepsSimulated+preTimeStepsSimulated);
            }
        }
    }

    for (unsigned int i=0; i<pDistrict->getnDECs(); i++) { // Cognet: Added this.
        bytes += pDistrict->getDEC(i)->nbFloatsRecorded()*sizeof(float)*(timeStepsSimulated+preTimeStepsSimulated);
    }
    return bytes;
}

size_t XmlScene::memoryUsageIrradiation() {

    // add the obstructing surfaces' 2 vectors of float (SW,LW)
    size_t bytes = pDistrict->getnSurfaces()*sizeof(float)*2*(timeStepsSimulated+preTimeStepsSimulated);
    // add the ground surfaces' 2 vectors of float (SW,LW)
    bytes += pDistrict->getnGrounds()*sizeof(float)*2*(timeStepsSimulated+preTimeStepsSimulated);
    // add the tree surfaces results
    for (size_t treeIndex=0; treeIndex<pDistrict->getnTrees(); ++treeIndex)
        bytes += pDistrict->getTree(treeIndex)->getnSurfaces()*sizeof(float)*1*(timeStepsSimulated+preTimeStepsSimulated); // SW
    // building results
    for (size_t buildingIndex=0; buildingIndex<pDistrict->getnBuildings(); ++buildingIndex) {
        for (size_t zoneIndex=0; zoneIndex<pDistrict->getBuilding(buildingIndex)->getnZones(); ++zoneIndex) {
            for (size_t wallIndex=0; wallIndex<pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getnWalls(); ++wallIndex) {
                // for each wall sum of 5 vectors of float (SW,DLout,3*DL)
                bytes += sizeof(*(pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getWall(wallIndex)));
                bytes += sizeof(float)*5*(timeStepsSimulated+preTimeStepsSimulated);
            }
            for (size_t roofIndex=0; roofIndex<pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getnRoofs(); ++roofIndex) {
                // for each wall sum of 2 vectors of float (SW,DLout)
                bytes += sizeof(*(pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getRoof(roofIndex)));
                bytes += sizeof(float)*2*(timeStepsSimulated+preTimeStepsSimulated);
            }
            for (size_t surfaceIndex=0; surfaceIndex<pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getnSurfaces(); ++surfaceIndex) {
                // for each wall sum of 2 vectors of float (SW,DLout)
                bytes += sizeof(*(pDistrict->getBuilding(buildingIndex)->getZone(zoneIndex)->getSurface(surfaceIndex)));
                bytes += sizeof(float)*2*(timeStepsSimulated+preTimeStepsSimulated);
            }
        }
    }
    return bytes;
}
