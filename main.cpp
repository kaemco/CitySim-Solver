/*
CitySim - http://leso.epfl.ch/citysim
Copyright (C) (2009-2013)  EPFL
*/

#include "scene.h"
#include "result.h"
#include "climate.h"

#ifdef GLUT
#include "GL/glut.h"
#endif

// JK - starting again the main of CitySim 03/02/09
// JK - removed all references to André's radiation model

using namespace std;

// *** main, CitySim         *** //
// *** jerome.kaempf@epfl.ch *** //

string about =
"\nCitySim Solver  Copyright (C) 2009-" + string(__DATE__).substr(11-4,4) + "  EPFL Solar Energy and Building Physics Laboratory (LESO-PB), kaemco LLC\n" +
"Build " + __DATE__ + " @ " + __TIME__ + "\n"
"This program comes with a LICENSE; for details type `CitySim --license'.\n\n";

string license = string("Redistribution and use ") +
"in source and binary forms, with or without modification, are permitted provided that the following conditions are met:\n\n" +
"1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.\n\n" +
"2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.\n\n" +
"3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.\n\n" +
"THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.";

int main(int argc, char *argv[])
{
    // print the disclaimer
    cout << about << endl;
    if (argc == 1) return 0;

    // graphical window (if needed)
#ifdef GLUT
    // load the glut
	glutInit(&argc, argv);
	glutInitWindowPosition( 0, 0 );
	glutInitWindowSize( 512u, 512u );
	glutInitDisplayMode( GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH ); //

	// Main window
	int win = glutCreateWindow("CitySim OpenGL Viewer");
#endif

    // starts the main
    ofstream logFileStream; // needed to have a log file
    ostream logStream(std::cerr.rdbuf()); // needed to choose between log file or standard output

    // JK - loading of the input file
    // first case: RAD file -> ViewFactorsCalculation, only irradiance on surfaces
    // second case: XML file -> AllCalculation

    try {

        string firstParameter;
        if (argc > 1) firstParameter = string(argv[1]);
        string secondParameter;
        if (argc > 2) secondParameter = string(argv[2]);
        string thirdParameter;
        if (argc > 3) thirdParameter = string(argv[3]);

        // prints the WARRANTY extract
        if (firstParameter == "--license") { cout << license << endl; return 0; }

        // go to the simulation bits
        if ( firstParameter == "-Ie" && secondParameter.size() > 4 && secondParameter.substr(secondParameter.size()-4,4) == ".rad" &&
             thirdParameter.size() > 4 && thirdParameter.substr(thirdParameter.size()-4,4) == ".cli" ) {

            cerr << "Radiance description file: " << firstParameter << endl;
            cerr << "Climate file: " << secondParameter << endl;
            cerr << "Comparison of the external irradiance." << endl;

            // loading the scene file
            Radscene radianceScene(firstParameter,secondParameter);

            // clear the results file
            remove( (firstParameter.substr(0,firstParameter.size()-4) + ".dat").c_str() );

            // writes the header in the output file
            radianceScene.writeHeader( firstParameter.substr(0,firstParameter.size()-4) + ".dat" );

            // loop on the first days
            for (unsigned int day=1;day<=1;day++) {
                for (unsigned int hour=9;hour<=17;hour++) {

                    // computes the shortwave using the SRA model
                    radianceScene.computeShortWave(day,hour);
                    // compare with Radiance by starting the executable
                    radianceScene.compareWithRadianceExternalIrradiance(day,hour);
                    radianceScene.writeResults( firstParameter.substr(0,firstParameter.size()-4) + ".dat" );
                    // save the total as an .out file
                    radianceScene.exportSWFile( firstParameter.substr(0,firstParameter.size()-4) + "_" + toString(day) + "_" + toString(hour) + "_SW.out" );
                    // clears the results
                    radianceScene.clearResults();

                }
            }
            return 0;
        }
        else if ( firstParameter.size() > 4 && firstParameter.substr(firstParameter.size()-4,4) == string(".rad") &&
             secondParameter.size() > 4 && secondParameter.substr(secondParameter.size()-4,4) == string(".rad") ) {

            cerr << "Radiance description file: " << firstParameter << endl;
            cerr << "Climate file: " << secondParameter << endl;

            // loading the scene file
            Radscene radianceScene(firstParameter,string(""));

            // exports the INP file
            radianceScene.exportInpFile();

            // do the calculation with Radiance for the given sky, creates the OUT file
            radianceScene.withRadianceClimate(secondParameter);

            return 0;
        }
        else if ( firstParameter.size() > 4 && firstParameter.substr(firstParameter.size()-4,4) == string(".rad") &&
                  secondParameter == string("dxf") ) {

                cerr << "Radiance description file: " << firstParameter << endl;
                cerr << "Converting to DXF..." << endl;

                // loading the scene file
                Radscene radianceScene(firstParameter,string(""));

                // export the .inp file
                radianceScene.exportDXF();

                return 0;
        }
        else if (   ((firstParameter.size() > 4) && (firstParameter.substr(firstParameter.size()-4,4) == string(".xml")) && (secondParameter.empty()))
                 || ((firstParameter == string("-q")) && (secondParameter.size() > 4) && (secondParameter.substr(secondParameter.size()-4,4) == string(".xml")))  ) {

            if (firstParameter == string("-q")) {
                // if -q then quiet output, redirected in files
                firstParameter = secondParameter;
                logFileStream.open(string(firstParameter.substr(0,firstParameter.size()-4) + ".log").c_str());
                if (logFileStream.fail())
                    throw(string("Unable to connect logFileStream to the .log file."));
                logStream.rdbuf(logFileStream.rdbuf());
                // As some logs using "cerr << " could not be directed in logStream, redirect cerr in the .log file too.
                std::cerr.rdbuf(logFileStream.rdbuf());
                //if (freopen(string(firstParameter.substr(0,firstParameter.size()-4) + ".log").c_str(), "a", stderr) == NULL)
                //    throw(string("Unable to redirect stderr to a file."));
            }

            // full simulation with the .xml file
            logStream << "XML description file: " << firstParameter << endl;

            // creates the XmlScene, which will populate under classes (pre-process)
            XmlScene xmlscene(firstParameter, &logStream);

            #ifdef DEBUG
            // outputs the radFile (for visualisation), the Inp file and the cumulative radiance files, plus a DXF file
            xmlscene.exportRadFile();
            xmlscene.exportRadFile(firstParameter.substr(0,firstParameter.size()-4) + "_triangulated.rad", true);
            xmlscene.exportInpFile();
            #endif
            //xmlscene.computeCumulativeRadiance();
            xmlscene.exportDXF();
            xmlscene.exportSTL();
            //cerr << "Writing XML file..." << endl;
            //xmlscene.exportXMLFile();

            // prepare the output files (headers)
            xmlscene.writeSWHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out");
            xmlscene.writeSWvHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_SWv.out");
            xmlscene.writeDLHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_DL.out");
            xmlscene.writeLWHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_LW.out");
            xmlscene.writeTHHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_TH.out");
            if (Model::thermalExplicit)
                xmlscene.writeTHExplicitHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_THExplicit.out");
            xmlscene.writeTSHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_TS.out");
            xmlscene.writeHCHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_HC.out");
            xmlscene.writeETHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_ET.out");
            xmlscene.writeCMHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_CM.out");

            // pre-simulation process
            xmlscene.computeViewFactors();
            xmlscene.writeVFText(firstParameter.substr(0,firstParameter.size()-4) + "_VF.out");

            // simulation process
            xmlscene.simulate();

            // save the results
            xmlscene.writeSWResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out");
            xmlscene.writeSWvResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_SWv.out");
            xmlscene.writeDLResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_DL.out");
            xmlscene.writeLWResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_LW.out");
            xmlscene.writeTHResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_TH.out");
            if (Model::thermalExplicit)
                xmlscene.writeTHExplicitResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_THExplicit.out");
            xmlscene.writeTSResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_TS.out");
            xmlscene.writeHCResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_HC.out");
            xmlscene.writeETResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_ET.out");
            xmlscene.writeCMResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_CM.out");

            // save climatic data
            xmlscene.writeClimaticDataText(firstParameter.substr(0,firstParameter.size()-4) + "_ClimaticData.out");
            xmlscene.writeAreaText(firstParameter.substr(0,firstParameter.size()-4) + "_Area.out");
            xmlscene.writeInertiaText(firstParameter.substr(0,firstParameter.size()-4) + "_Inertia.out");

            // save aggregated results
//            xmlscene.writeMonthlyResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_MonthlyResults.out");
            xmlscene.writeYearlyResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_YearlyResults.out");
            xmlscene.writeYearlyResultsPerBuildingText(firstParameter.substr(0,firstParameter.size()-4) + "_YearlyResultsPerBuilding.out");

            // export in GML format
            xmlscene.exportGML(firstParameter.substr(0,firstParameter.size()-4) + ".gml");

            return 0;

        }
        else if ( firstParameter.size() > 1 && firstParameter == string("-I") && !secondParameter.empty() ) {

            // only irradiation calculation with the .xml file, the secondParameter contains the firstParameter
            firstParameter = secondParameter;
            cerr << "XML description file: " << firstParameter << endl;

            // creates the XMLscene, which will populate under classes (pre-process)
            XmlScene xmlscene(firstParameter);

            // outputs the radFile (for visualisation)
            xmlscene.exportRadFile();
            xmlscene.exportRadFile(firstParameter.substr(0,firstParameter.size()-4) + "_triangulated.rad");
            xmlscene.exportInpFile(firstParameter.substr(0,firstParameter.size()-4) + ".inp", true);

            // prepare the output files (headers)
            xmlscene.writeSWHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out");

            // pre-simulation process
            xmlscene.computeViewFactors();
            xmlscene.writeVFText(firstParameter.substr(0,firstParameter.size()-4) + "_VF.out");

            // simulation process
            xmlscene.simulateRadiation(false); // do not do the daylight calculation

            // save the results
            xmlscene.writeSWResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out");

            return 0;

        }
        else if ( firstParameter.size() > 1 && firstParameter == string("-Id") && !secondParameter.empty() ) {

            // only irradiation calculation with the .xml file, the secondParameter contains the firstParameter
            firstParameter = secondParameter;
            cerr << "XML description file: " << firstParameter << endl;

            // creates the XMLscene, which will populate under classes (pre-process)
            XmlScene xmlscene(firstParameter);

            // outputs the radFile (for visualisation)
            xmlscene.exportRadFile();
            xmlscene.exportRadFile(firstParameter.substr(0,firstParameter.size()-4) + "_triangulated.rad");
            xmlscene.exportInpFile();

            // prepare the output files (headers)
            xmlscene.writeSWHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out");
            xmlscene.writeSWvHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_SWv.out");
            xmlscene.writeDLHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_DL.out");

            // pre-simulation process
            xmlscene.computeViewFactors();
            xmlscene.writeVFText(firstParameter.substr(0,firstParameter.size()-4) + "_VF.out");

            // simulation process
            xmlscene.simulateRadiation();

            // save the results
            xmlscene.writeSWResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out");
            xmlscene.writeSWvResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_SWv.out");
            xmlscene.writeDLResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_DL.out");

            return 0;

        }
        else if ( firstParameter.size() > 1 && firstParameter == string("-Ia")
                  && (secondParameter.size() > 4) && (secondParameter.substr(secondParameter.size()-4,4) == string(".xml")) ) {

            // only irradiation calculation with the .xml file, the secondParameter contains the firstParameter
            firstParameter = secondParameter;
            cerr << "XML description file: " << firstParameter << endl;

            // creates the XMLscene, which will populate under classes (pre-process)
            XmlScene xmlscene(firstParameter);

            // outputs the radFile (for visualisation)
            xmlscene.exportRadFile();
            xmlscene.exportRadFile(firstParameter.substr(0,firstParameter.size()-4) + "_triangulated.rad");
            xmlscene.exportInpFile();

            // initialise the far field
            xmlscene.initialiseFarField();

            // prepare the annual climate sky
            cout << "Exporting annual sky and ground." << endl;
            xmlscene.exportCumulativeRadiance();

            // prepare the output files (headers)
            xmlscene.writeSWHeaderText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out", "Irradiation(Wh/m2)");

            // pre-simulation process
            xmlscene.computeViewFactors();
            xmlscene.writeVFText(firstParameter.substr(0,firstParameter.size()-4) + "_VF.out");

            // simulation process from:
            xmlscene.simulateCumulativeIrradiance();

            // save the results
            xmlscene.writeSWResultsText(firstParameter.substr(0,firstParameter.size()-4) + "_SW.out");

            return 0;

        }
        else if ( firstParameter.size() > 1 && firstParameter == string("-Ih")
                  && (secondParameter.size() > 4) && (secondParameter.substr(secondParameter.size()-4,4) == string(".xml")) ) {

            // only irradiation calculation with the .xml file, the secondParameter contains the firstParameter
            firstParameter = secondParameter;
            cerr << "XML description file: " << firstParameter << endl;

            // creates the XMLscene, which will populate under classes (pre-process)
            XmlScene xmlscene(firstParameter);

            // outputs the radFile (for visualisation)
            xmlscene.exportRadFile();
            xmlscene.exportRadFile(firstParameter.substr(0,firstParameter.size()-4) + "_triangulated.rad");
            xmlscene.exportInpFile();

            // initialise the far field
            xmlscene.initialiseFarField();

            // prepare the annual climate sky
            cout << "Exporting hourly sky and ground." << endl;
            xmlscene.exportHourlyRadiance();

            return 0;

        }
        else if ( firstParameter.size() > 4 && firstParameter.substr(firstParameter.size()-4,4) == string(".xml") &&
                  secondParameter == string("dxf") ) {

            cerr << "XML description file: " << firstParameter << endl;
            cerr << "Converting to DXF..." << endl;

            // creates the XmlScene, which will populate under classes (pre-process)
            XmlScene xmlscene(firstParameter);

            // exports the DXF file
            xmlscene.exportDXF();

            return 0;
        }
        else if ( firstParameter.size() > 4 && firstParameter.substr(firstParameter.size()-4,4) == string(".cli") &&
                  secondParameter.empty() ) {

            cerr << "Climate file: " << firstParameter << endl;

            Scene dummyScene;
            // creates an empty scene with the climate, and output the maxwell values if G_h is given with the debug version
            Scene scene(firstParameter, firstParameter); // defines the climate as second parameter, and the inputfile as the first parameter
            // creates the cumulative radiance file
            scene.exportCumulativeRadiance();

            return 0;
        }
        else if (firstParameter == string("-M") && (secondParameter.size() > 4) && (secondParameter.substr(secondParameter.size()-4,4) == string(".xml"))) {

            // full simulation with the .xml file
            logStream << "XML description file: " << secondParameter << endl;

            // creates the XmlScene, which will populate under classes (pre-process)
            XmlScene xmlscene(secondParameter);

            xmlscene.exportDXF();
            xmlscene.exportSTL();

            // pre-simulation process
            xmlscene.computeViewFactors();

            // open text file for the Matrix of View Factors
            fstream textFile((secondParameter.substr(0,secondParameter.size()-4) + "_MA.out").c_str(),ios::out|ios::binary|ios::trunc);
            // write the matrix itself
            for (unsigned int i=0; i < xmlscene.getnAi(); ++i) { // loop on the number of elements
                for (unsigned int index=xmlscene.getAi(i); index < xmlscene.getAi(i+1); ++index) {
                    textFile << i << "\t" << xmlscene.getAj(index) << "\t" << xmlscene.getAn(index)/(((Surface*)(xmlscene.getDATARadiationScene()->GetSurface(xmlscene.getAj(index)).SurfaceDelegate()))->getShortWaveReflectance())*M_PI << endl;
                }
            }
            textFile.close();

            // open text file for the Matrix of distances
            textFile.open((secondParameter.substr(0,secondParameter.size()-4) + "_Md.out").c_str(),ios::out|ios::binary|ios::trunc);
            // write the matrix itself
            for (unsigned int surfaceIndex_1=0; surfaceIndex_1<xmlscene.getDATARadiationScene()->SurfaceCount(); ++surfaceIndex_1) {
                for (unsigned int surfaceIndex_2=0; surfaceIndex_2<xmlscene.getDATARadiationScene()->SurfaceCount(); ++surfaceIndex_2) {
                    GENPoint distance = ((Surface*)(xmlscene.getDATARadiationScene()->GetSurface(surfaceIndex_1).SurfaceDelegate()))->computeCentroid()
                                        -((Surface*)(xmlscene.getDATARadiationScene()->GetSurface(surfaceIndex_2).SurfaceDelegate()))->computeCentroid();
                    textFile << surfaceIndex_1 << "\t" << surfaceIndex_2 << "\t" << sqrt(dot_product(distance,distance)) << endl;
                }
            }
            textFile.close();

            return 0;
        }
        else throw(string("Unknown file format: ") + firstParameter);

    }
    catch(string msg) { logStream << "\n(Caught) " << msg << endl; return 1; }
    catch(exception &e) { logStream << "\n(Caught Standard Exception) " << e.what() << endl; return 1; }
    catch(...) { logStream << "\nUnhandled Exception" << endl; return 1; }

    return 0;
}
