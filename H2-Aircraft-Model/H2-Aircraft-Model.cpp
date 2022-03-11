// H2-Aircraft-Model.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "headers/AircraftModel.h"
#include "headers/MathTools.h"
#include "headers/PlaneMakerTools.h"

int main()
{
    // ACF File Location
    std::string acf_filepath =
        "D:\\Games\\X-Plane 11\\Aircraft\\Extra Aircraft\\ATR72-500\\ATR72.acf";

    PlaneMakerTools::set_engine_data(0, 0, 0, acf_filepath);

    double y = AircraftModel::compute_cruise_BSFC_PW127(5);
    y = y * 3.6e+9;
    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
