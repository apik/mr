#ifndef _LOGGER_HPP_
#define _LOGGER_HPP_

#include <iostream>
#include <sstream>

/* consider adding boost thread id since we'll want to know whose
   writting and
   * won't want to repeat it for every single call */

/* consider adding policy class to allow users to redirect logging to
   specific
   * files via the command line
   */

enum loglevel_e
  {logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4};

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
  logIt & operator<<(T const & value)
  {
    _buffer << value;
    return *this;
  }
  logIt & operator<<(std::ostream& (*fun)(std::ostream&)) { std::cout << std::endl; return *this;}
  
  ~logIt()
  {
    _buffer << std::endl;
    // This is atomic according to the POSIX standard
    // http://www.gnu.org/s/libc/manual/html_node/Streams-and-Threads.html
    std::cerr << _buffer.str();
  }
  
private:
  std::ostringstream _buffer;
};

extern loglevel_e loglevel;

#define lout(level)       \
  if (level > loglevel) ; \
  else logIt(level)

#endif


