//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef _LOGGER_HPP_
#define _LOGGER_HPP_

#include <iostream>
#include <sstream>

namespace mr
{

  enum loglevel_e
    {logERROR, logWARNING, logINFO, logDEBUG};
  
  class logIt
  {
  public:
    logIt(loglevel_e _loglevel = logERROR) {

      if (_loglevel == logERROR)
        _buffer << "[ERROR]:    " ;
      else if (_loglevel == logWARNING)
        _buffer << "[WARNING]:  " ;
      else if (_loglevel == logINFO)
        _buffer << "[INFO]:     " ;
      else if (_loglevel == logDEBUG)
        _buffer << "[DEBUG]:    " ;
      else
        {
          _buffer << "[DBG" <<  (_loglevel - logDEBUG) << "]"
                  << std::string(
                                 _loglevel > logDEBUG
                                 ? (_loglevel - logDEBUG) * 4
                                 : 1
                                 , ' ');
        }
    }
  
    template <typename T>
    logIt & operator << (T const & value)
    {
      _buffer << value;
      return *this;
    }
    logIt & operator << (std::ostream& (*fun)(std::ostream&)) { std::cout << std::endl; return *this;}
  
    ~logIt()
    {
      _buffer << std::endl;
      std::cerr << _buffer.str();
    }
  
  private:
    std::ostringstream _buffer;
  };

  extern loglevel_e loglevel;
  
} // namespace mr



#define lout(level)                             \
  if (level > loglevel) ;                       \
  else logIt(level)
#endif


