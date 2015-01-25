/*******************************************************************************
 * Copyright (c) 2015 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the GNU LGPL v3
 *******************************************************************************
 *
 * Filename: main.cpp
 *
 * Description:
 *      description
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2015-01-15   wm              Initial version
 *
 ******************************************************************************/

#include "ABCSpeedup.hpp"

#include <fstream>
#include <string>
#include <iostream>
#include <cassert>
#include <cstddef>
#include <algorithm>
#include <utility>

double Accuracy(std::vector<std::string> && refvs, std::vector<std::string> && testvs)
{
    assert(refvs.size() == testvs.size());

    const std::size_t nelem = refvs.size();

    std::vector<std::size_t> refdata(nelem);
    std::vector<std::size_t> testdata(nelem);

    std::transform(refvs.cbegin(), refvs.cend(), refdata.begin(),
        [](const std::string & str) -> std::size_t
        {
            return atoi(str.c_str());
        }
    );

    std::transform(testvs.cbegin(), testvs.cend(), testdata.begin(),
        [](const std::string & str) -> std::size_t
        {
            return atoi(str.c_str());
        }
    );

    std::size_t correct = 0;
    std::size_t baseCorrect = 0;

    const std::size_t * refdata_p = refdata.data();
    const std::size_t * testdata_p = testdata.data();

    for (std::size_t i = 0; i < nelem; ++i)
    {
        for (std::size_t j = i + 1; j < nelem; ++j)
        {
            if ((refdata_p[i] == refdata_p[j]) == (testdata_p[i] == testdata_p[j]))
            {
                correct += 1;
            }
            if ((refdata_p[i] == refdata_p[j]) == (i == j))
            {
                baseCorrect += 1;
            }
        }
    }

    double accuracy = (double)std::max(correct - baseCorrect, 0u) / (nelem * (nelem - 1) / 2 - baseCorrect);

    return accuracy;
}

int main(int argc, char **argv)
{
    // main input validation
    if (argc == 3)
    {
        std::ifstream fjson(argv[1]);
        std::vector<std::string> vjson;

        for (std::string line; std::getline(fjson, line);)
        {
            vjson.push_back(line);
        }

        ABCSpeedup worker;
        std::vector<std::string> assignments = worker.cluster(vjson);

        std::ifstream frefdata(argv[2]);
        std::vector<std::string> vrefdata;

        for (std::string line; std::getline(frefdata, line);)
        {
            vrefdata.push_back(line);
        }

        std::cerr << "Ref # " << vrefdata.size() << ", Got # " << assignments.size() << std::endl;
        assert(vrefdata.size() == assignments.size());

        double accuracy = Accuracy(std::move(vrefdata), std::move(assignments));

        std::cout << "Accuracy: [" << accuracy << "]" << std::endl;
    }
    else
    {

    }

    return 0;
}
