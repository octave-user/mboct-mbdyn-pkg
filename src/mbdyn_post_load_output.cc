//  Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; If not, see <http://www.gnu.org/licenses/>.

#include "config.h"

#define NDEBUG
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#ifndef NDEBUG
#define BOUNDS_CHECKING
#endif

#include <octave/oct.h>

// PKG_ADD: autoload ("mbdyn_post_load_output_out", "__mboct_mbdyn__.oct");
// PKG_DEL: autoload ("mbdyn_post_load_output_out", "__mboct_mbdyn__.oct", "remove");

DEFUN_DLD(mbdyn_post_load_output_out, args, nargout,
          "-*- texinfo -*-\n"
          "@deftypefn {} [@var{Time}, @var{TStep}, @var{NIter}, @var{ResErr}, @var{SolErr}, @var{SolConv}, @var{OutputFlag}] = mbdyn_post_load_output_out(@var{mbdyn_filename}, @var{append_rows}, @var{auto_resize_rows})\n\n"
          "Load MBDyn's .out file from disk.\n\n"
          "@var{Time} @dots{} Simulation time.\n\n"
          "@var{TStep} @dots{} Time step.\n\n"
          "@var{NIter} @dots{} Number of iterations.\n\n"
          "@var{ResErr} @dots{} Residual error.\n\n"
          "@var{SolErr} @dots{} Solution error.\n\n"
          "@var{SolConv} @dots{} Convergence on solution.\n\n"
          "@var{OutputFlag} @dots{} A one indicates, that this time step has been written to all output files (e.g. defined by output meter).\n\n"
          "@var{mbdyn_filename} @dots{} MBDyn output filename excluding extension (e.g. defined by -o flag).\n\n"
          "@var{append_rows} @dots{} Hint for memory reallocation.\n\n"
          "@var{auto_resize_rows} @dots{} Automatically resize the output (e.g. if the file is read while MBDyn is still writing to it).\n\n"
          "@end deftypefn")
{
    octave_value_list retval;

    octave_idx_type nargin = args.length();

    if (nargin < 1 || nargin > 3)
    {
        print_usage();
        return retval;
    }

    if (!args(0).is_string() && !args(0).is_sq_string())
    {
        error_with_id("mbdyn:post", "mbdyn_mov_filename must be a string!");
        return retval;
    }

    const std::string out_filename = args(0).string_value() + ".out";

    octave_idx_type append_rows = 1024;

    if (nargin >= 2)
    {
        if (!args(1).is_scalar_type())
        {
            error_with_id("mbdyn:post", "append_rows must be a scalar value!");
            return retval;
        }

        append_rows = args(1).int_value();
    }

    if (append_rows <= 1)
    {
        error_with_id("mbdyn:post", "append_rows <= 1");
        return retval;
    }

    bool auto_resize_rows = true;

    if (nargin >= 3)
    {
        if (!args(2).OV_ISLOGICAL() || !args(2).is_scalar_type())
        {
            error_with_id("mbdyn:post", "auto_resize_rows must be a boolean scalar");
            return retval;
        }

        auto_resize_rows = args(2).bool_value();
    }

    const bool fTStep = nargout >= 2;
    const bool fNIter = nargout >= 3;
    const bool fResErr = nargout >= 4;
    const bool fSolErr = nargout >= 5;
    const bool fSolConv = nargout >= 6;
    const bool fOutputFlag = nargout >= 7;

    std::ifstream out_file;

    out_file.open(out_filename.c_str());

    if (!out_file.good())
    {
        error_with_id("mbdyn:post", "could not open file \"%s\"!", out_filename.c_str());
        return retval;
    }

    octave_idx_type N = append_rows;

    out_file.seekg(0, std::ios::end);

    std::ios::off_type pos = out_file.tellg();
    bool bGotNumLines = false;
    std::string sStep;

    for (std::ios::off_type i = pos - 1; i >= 0; --i)
    {
        char c = '\0';

        out_file.seekg(i, std::ios::beg);

        if (!out_file.get(c))
        {
            error_with_id("mbdyn:post", "failed to read from file \"%s\"!", out_filename.c_str());
            return retval;
        }

        if (c == '\n')
        {
            octave_idx_type iStep;

            out_file >> sStep >> iStep;

            if (!out_file)
            {
                out_file.clear();
                continue;
            }

            if (sStep == "Step")
            {
                N = iStep + 1;
                bGotNumLines = true;
                break;
            }
        }
    }

    if (!bGotNumLines)
    {
        error_with_id("mbdyn:post", "failed to extract the number of lines from file \"%s\"!", out_filename.c_str());
        return retval;
    }

    out_file.seekg(0, std::ios::beg);

    ColumnVector Time(N, 0.);
    ColumnVector TStep(fTStep ? N : 0, 0.);
    int32NDArray NIter(dim_vector(fNIter ? N : 0, 1), 0);
    ColumnVector ResErr(fResErr ? N : 0, 0.);
    ColumnVector SolErr(fSolErr ? N : 0, 0.);
    int32NDArray SolConv(dim_vector(fSolConv ? N : 0, 1), 0);
    boolNDArray OutputFlag(dim_vector(fOutputFlag ? N : 0, 1), false);

    std::string line, strStep;
    std::istringstream is;

    octave_idx_type i = -1, iLine = 0;
    octave_idx_type iPrevStep = -1;

    while (true)
    {
        OCTAVE_QUIT;

        getline(out_file, line, '\n');

        if (out_file.eof())
        {
            break;
        }
        else if (!out_file.good())
        {
            error_with_id("mbdyn:post", "during reading file \"%s\"!", out_filename.c_str());
            return retval;
        }

        ++iLine;
        is.clear();
        is.str(line);

        // # Step Time TStep NIter ResErr SolErr SolConv
        // Step 3094 0.343636 6.1704e-13 3 6.07262e-07 0 0 0
        octave_idx_type Step;

        is >> strStep >> Step;

        if (strStep != "Step")
        {
            // If we perform an eigenanalysis the eigenvalues are dumped in the out file
            // Therefore continue if the line does not begin with "Step"
            continue;
        }

        if (Step != iPrevStep + 1)
        {
            // Do not return dummy steps!
            warning_with_id("mbdyn:post", "Step %Ld in file %s, line %Ld is out of sequence (previous step %Ld)", static_cast<long long>(Step), out_filename.c_str(), static_cast<long long>(iLine), static_cast<long long>(iPrevStep));
            continue;
        }

        if (++i >= N)
        {
            warning_with_id("mbdyn:post", "the file \"%s\" has been changed during last read operation!", out_filename.c_str());

            if (auto_resize_rows)
            {
                N += append_rows;

                Time.resize(N);

                if (fTStep)
                    TStep.resize(N);

                if (fNIter)
                    NIter.resize1(N);

                if (fResErr)
                    ResErr.resize(N);

                if (fSolErr)
                    SolErr.resize(N);

                if (fSolConv)
                    SolConv.resize1(N);

                if (fOutputFlag)
                    OutputFlag.resize1(N);
            }
            else
            {
                --i;
                break;
            }
        }

        is >> Time(i);

        if (fTStep)
            is >> TStep(i);

        if (fNIter)
            is >> NIter(i);

        if (fResErr)
            is >> ResErr(i);

        if (fSolErr)
            is >> SolErr(i);

        if (fSolConv)
            is >> SolConv(i);

        if (fOutputFlag)
            is >> OutputFlag(i);

        if (is.fail())
        {
            error_with_id("mbdyn:post", "during reading file \"%s\"\n\tline: \"%s\"", out_filename.c_str(), line.c_str());
            return retval;
        }

        iPrevStep = Step;
    }

    assert(i < N);

    N = i + 1;

    Time.resize(N);

    if (fTStep)
        TStep.resize(N);

    if (fNIter)
        NIter.resize1(N);

    if (fResErr)
        ResErr.resize(N);

    if (fSolErr)
        SolErr.resize(N);

    if (fSolConv)
        SolConv.resize1(N);

    if (fOutputFlag)
        OutputFlag.resize1(N);

    retval.append(octave_value(Time));

    if (fTStep)
        retval.append(octave_value(TStep));

    if (fNIter)
        retval.append(octave_value(NIter));

    if (fResErr)
        retval.append(octave_value(ResErr));

    if (fSolErr)
        retval.append(octave_value(SolErr));

    if (fSolConv)
        retval.append(octave_value(SolConv));

    if (fOutputFlag)
        retval.append(octave_value(OutputFlag));

    return retval;
}

class mbdyn_file_output
{
    public:
    inline mbdyn_file_output();
    inline void resize(octave_idx_type row_size_request, octave_idx_type column_size_request, octave_idx_type append_rows, octave_idx_type append_columns);
    octave_idx_type rows()const{ return current_row_size; }
    octave_idx_type columns()const{ return current_column_size; }
    double& operator()(octave_idx_type i, octave_idx_type j);
    const Matrix& get_data();
    private:
    octave_idx_type current_row_size;
    octave_idx_type current_column_size;
    Matrix data;
};

mbdyn_file_output::mbdyn_file_output()
    :current_row_size(0), current_column_size(0)
{

}

void mbdyn_file_output::resize(octave_idx_type row_size_request, octave_idx_type column_size_request, octave_idx_type append_rows, octave_idx_type append_columns)
{
    octave_idx_type capacity_rows = data.rows();
    octave_idx_type capacity_columns = data.columns();

    if (capacity_rows < row_size_request || capacity_columns < column_size_request)
    {
        while (capacity_rows < row_size_request)
            capacity_rows += append_rows;

        while (capacity_columns < column_size_request)
            capacity_columns += append_columns;

        data.resize(capacity_rows,capacity_columns);
    }

    current_row_size = row_size_request;
    current_column_size = column_size_request;

    assert(current_row_size <= data.rows());
    assert(current_column_size <= data.columns());
}

double& mbdyn_file_output::operator()(octave_idx_type i, octave_idx_type j)
{
    assert(i >= 0 && i < current_row_size);
    assert(j >= 0 && j < current_column_size);
    return data(i,j);
}

const Matrix& mbdyn_file_output::get_data()
{
    if (data.rows() > current_row_size || data.columns() > current_column_size)
    {
        data.resize(current_row_size, current_column_size);
    }

    return data;
}

// PKG_ADD: autoload ("mbdyn_post_load_output", "__mboct_mbdyn__.oct");
// PKG_DEL: autoload ("mbdyn_post_load_output", "__mboct_mbdyn__.oct", "remove");

DEFUN_DLD(mbdyn_post_load_output, args, nargout,
          "-*- texinfo -*-\n"
          "@deftypefn {} [@var{id}, @var{data}] = mbdyn_post_load_output(@var{filename}, @var{default_column_count}, @var{filter_id}, @var{append_rows}, @var{append_columns}, @var{output_step}, @var{auto_resize_rows})\n"
          "Load a generic output file from MBDyn.\n\n"
          "@var{filename} @dots{} Output file generated by MBDyn including extension.\n\n"
          "@var{column count} @dots{} Number of columns to load excluding the id column.\n\n"
          "@var{filter_id} @dots{} Load only those entities listed in this array. If empty, everything is loaded.\n\n"
          "@var{append rows} @dots{} Hint for memory reallocation.\n\n"
          "@end deftypefn")
{
    octave_value_list retval;

    const octave_idx_type nargin = args.length();

    if (nargin < 2 || nargin > 7)
    {
        print_usage();
        return retval;
    }

    if (!args(0).is_string() && !args(0).is_sq_string())
    {
        error_with_id("mbdyn:post", "filename must be a string!");
        return retval;
    }

    const std::string filename = args(0).string_value();

    if (!args(1).is_scalar_type())
    {
        error_with_id("mbdyn:post", "column_count must be a scalar value!");
        return retval;
    }

    const octave_idx_type default_column_count = args(1).int_value();

    if (default_column_count < 0)
    {
        error_with_id("mbdyn:post", "column_count must be greater than zero!");
        return retval;
    }

    int32NDArray filter_id;

    if (nargin >= 3)
    {
        if (!args(2).is_matrix_type() && !args(2).is_scalar_type() && !args(2).is_range())
        {
            error_with_id("mbdyn:post", "filter_id must be a matrix type!");
            return retval;
        }

        if ((args(2).rows() == 1 && args(2).columns() >= 1) ||
            (args(2).rows() >= 1 && args(2).columns() == 1))
        {
            filter_id = args(2).int32_array_value();
        }
        else if (args(2).rows() > 1 && args(2).columns() > 1)
        {
            error_with_id("mbdyn:post", "filter_id must be a vector!");
            return retval;
        }
    }

    std::set<int32_t> filter_id_set;

    for (octave_idx_type i = 0; i < filter_id.numel(); ++i)
        filter_id_set.insert(filter_id(i));

    octave_idx_type append_rows = 1024;

    if (nargin >= 4)
    {
        if (!args(3).is_scalar_type())
        {
            error_with_id("mbdyn:post", "append_rows must be scalar!");
            return retval;
        }

        append_rows = args(3).int_value();

        if (append_rows <= 0)
        {
            error_with_id("mbdyn:post", "append_rows must be greater than zero!");
            return retval;
        }
    }

    octave_idx_type append_columns = 6;

    if (nargin >= 5)
    {
        if (!args(4).is_scalar_type())
        {
            error_with_id("mbdyn:post", "append_columns must be scalar!");
            return retval;
        }

        append_columns = args(4).int_value();

        if (append_columns <= 0)
        {
            error_with_id("mbdyn:post", "append_columns must be greater than zero!");
            return retval;
        }
    }

    octave_idx_type output_step = 1;

    if (nargin >= 6)
    {
        if (!args(5).is_scalar_type())
        {
            error_with_id("mbdyn:post", "output_step must be a scalar!");
            return retval;
        }

        output_step = args(5).int_value();

        if (output_step < 1)
        {
            error_with_id("mbdyn:post", "output step must be greater than zero!");
        }
    }

    bool auto_resize_rows = true;

    if (nargin >= 7)
    {

        if (!args(6).OV_ISLOGICAL() || !args(6).is_scalar_type())
        {
            error_with_id("mbdyn:post", "auto_resize_rows must be a boolean scalar value");
            return retval;
        }

        auto_resize_rows = args(6).bool_value();
    }

    std::ifstream mbdyn_file;

    mbdyn_file.open(filename.c_str());

    if (!mbdyn_file.good())
    {
        error_with_id("mbdyn:post", "could not open file \"%s\"!", filename.c_str());
        return retval;
    }

    typedef std::map<int32_t, mbdyn_file_output> mbdyn_map_t;

    typedef mbdyn_map_t::iterator mbdyn_iter_t;

    mbdyn_map_t output_data;

    std::string line;
    std::istringstream is;
    octave_idx_type inum_elements = 0;
    octave_idx_type ielement = 0;
    octave_idx_type itime_step = 0;

    while (true)
    {
        OCTAVE_QUIT;

        std::getline(mbdyn_file, line,'\n');

        if (mbdyn_file.eof())
        {
            break;
        }
        else if (!mbdyn_file.good())
        {
            error_with_id("mbdyn:post", "while reading file \"%s\"!", filename.c_str());
            return retval;
        }

        if (line.length() > 0 && line[line.length() - 1] == '\r') // If we are running octave on Cygwin and mbdyn on mingw32 we have to consider the \r\n line endings
        {
            line.resize(line.length() - 1);
        }

        is.clear();
        is.str(line);

        int32_t current_id = -1;

        is >> current_id;

        if (!is.good())
        {
            error_with_id("mbdyn:post", "could not read id number from file \"%s\"!", filename.c_str());
            return retval;
        }

        bool foutput_current_id = false;

        if (filter_id_set.size() != 0)
        {
            foutput_current_id = filter_id_set.find(current_id) != filter_id_set.end();
        }
        else
        {
            // if no node filter is specified output every node!
            foutput_current_id = true;
        }

        if (foutput_current_id)
        {
            if (itime_step == 0)
            {
                if (output_data.end() == output_data.find(current_id))
                {
                    ++inum_elements;
                }
                else
                {
                    ++itime_step;
                }
            }
            else
            {
                if (++ielement > inum_elements)
                {
                    ielement = 1;
                    ++itime_step;
                }
            }

            if (itime_step % output_step)
            {
                continue;
            }

            mbdyn_file_output& node_output = output_data[current_id];

            octave_idx_type current_row_index = node_output.rows();

            for (octave_idx_type current_column_index = 0; is.good(); ++current_column_index)
            {
                double val;

                is >> val;

                if (is.fail())
                {
                    // replace inf and nan with NAN
                    val = NAN;
                    is.clear();
                    std::string token;
                    is >> token;
                }

                const octave_idx_type row_count = std::max(node_output.rows(), current_row_index + 1);
                const octave_idx_type column_count = std::max(node_output.columns(), current_column_index + 1);

                if (!auto_resize_rows && current_row_index >= append_rows)
                {
                    goto exit_load_data_loop;
                }

                node_output.resize(row_count, column_count, append_rows, append_columns);
                node_output(current_row_index, current_column_index) = val;
            }
        }
    }

exit_load_data_loop:
    int32NDArray id;
    id.resize1(nargout >= 1 ? output_data.size() : 0);
    Cell data(nargout >= 2 ? output_data.size() : 0, 1, Matrix());

    octave_idx_type i = 0;

    for (mbdyn_iter_t it = output_data.begin(); it != output_data.end(); ++it, ++i)
    {
        if (nargout >= 1)
            id(i) = it->first;

        if (nargout >= 2)
            data(i) = octave_value(it->second.get_data());
    }

    if (nargout >= 1)
        retval.append(octave_value(id));

    if (nargout >= 2)
        retval.append(octave_value(data));

    return retval;
}
