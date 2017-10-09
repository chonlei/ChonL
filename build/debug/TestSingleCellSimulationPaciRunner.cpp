/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "projects/ChonL/test/TestSingleCellSimulationPaci.hpp"

static TestSingleCellParameterSweep suite_TestSingleCellParameterSweep;

static CxxTest::List Tests_TestSingleCellParameterSweep = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestSingleCellParameterSweep( "projects/ChonL/test/TestSingleCellSimulationPaci.hpp", 89, "TestSingleCellParameterSweep", suite_TestSingleCellParameterSweep, Tests_TestSingleCellParameterSweep );

static class TestDescription_TestSingleCellParameterSweep_TestPaciSimulation : public CxxTest::RealTestDescription {
public:
 TestDescription_TestSingleCellParameterSweep_TestPaciSimulation() : CxxTest::RealTestDescription( Tests_TestSingleCellParameterSweep, suiteDescription_TestSingleCellParameterSweep, 92, "TestPaciSimulation" ) {}
 void runTest() { suite_TestSingleCellParameterSweep.TestPaciSimulation(); }
} testDescription_TestSingleCellParameterSweep_TestPaciSimulation;

#include <cxxtest/Root.cpp>
