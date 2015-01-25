/*******************************************************************************
 * Copyright (c) 2015 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the GNU LGPL v3
 *******************************************************************************
 *
 * Filename: timestamp.hpp
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
 * 2015-01-16   wm              Initial version
 *
 ******************************************************************************/

#ifndef TIMESTAMP_HPP_
#define TIMESTAMP_HPP_

#include <chrono>
#include <string>

struct Timestamp
{
    typedef std::chrono::steady_clock::time_point timepoint_type;
    explicit Timestamp(timepoint_type start)
    :
        m_start(start)
    {

    }

    std::string now(timepoint_type moment)
    {
//        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(moment - m_start).count();
//        return "[" + std::to_string(duration) + " us]";

        std::size_t duration = std::chrono::duration_cast<std::chrono::milliseconds>(moment - m_start).count();
        return "[" + std::to_string(duration) + " ms]";
    }

private:
    const timepoint_type m_start;
};

#endif /* TIMESTAMP_HPP_ */
