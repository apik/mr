#ifndef __GNUPLOT_HPP__
#define __GNUPLOT_HPP__
#include <string>
#include <fstream>

// Single graph
class Plot1
{
  std::string fname;
  std::ofstream  of;
  static const size_t fw = 20;         // Field width

public:
  // fname - w/o extension
  Plot1(const std::string& fname_, const std::string& title, const std::string& xlbl, const std::string& ylbl, const std::string& legend1) : fname(fname_)
  {
    // Plot script
    std::ofstream plt((fname + ".plt").c_str(), std::ofstream::out);

    plt << "# Run with gnuplot " << fname << ".plt" << std::endl;
    // plt << "set format y \"%.1f\"" << std::endl;
    plt << "set xtics  norangelimit" << std::endl;
    plt << "set grid x y" << std::endl;
    plt << "set mxtics" << std::endl;
    plt << "set mytics" << std::endl;
    plt << "set title \""<< title << "\"" << std::endl;
    plt << "set xlabel \"" << xlbl << "\"" << std::endl;
    plt << "set ylabel \"" << ylbl << "\"" << std::endl;
    plt << "label1 = \"" << legend1 << "\"" << std::endl;
    plt << "plot '" << fname << ".dat' using 1:2 with lines title label1 lt 1" << std::endl;
    plt << "" << std::endl;

    // Data file
    of.open((fname + ".dat").c_str());

    // Header
    of << "#" << std::setw(fw-3) << "X  " << std::setw(fw-3) <<  "Y1" << std::endl;
    of << std::scientific << std::setprecision(fw-10);
  }

  void add(long double x, long double y)
  {
    of << std::setw(fw) << x << std::setw(fw) << y << std::endl;
  }

}; 

class Plot2
{
  std::string fname;
  std::ofstream  of;
  static const size_t fw = 20;         // Field width

public:
  // fname - w/o extension
  Plot2(const std::string& fname_, const std::string& title, const std::string& xlbl, const std::string& ylbl, 
        const std::string& legend1, const std::string& legend2) : fname(fname_)
  {
    // Plot script
    std::ofstream plt((fname + ".plt").c_str(), std::ofstream::out);

    plt << "# Run with gnuplot " << fname << ".plt" << std::endl;
    plt << "set multiplot" << std::endl;
 
    // plt << "set format y \"%.1f\"" << std::endl;
    plt << "set xtics  norangelimit" << std::endl;
    plt << "set grid x y" << std::endl;
    plt << "set mxtics" << std::endl;
    plt << "set mytics" << std::endl;
    plt << "set title \""<< title << "\"" << std::endl;
    plt << "set xlabel \"" << xlbl << "\"" << std::endl;
    plt << "set ylabel \"" << ylbl << "\"" << std::endl;
    plt << "label1 = \"" << legend1 << "\"" << std::endl;
    plt << "label2 = \"" << legend2 << "\"" << std::endl;
    plt << "plot '" << fname << ".dat' using 1:2 with lines title label1 lt 1," 
        << "     '" << fname << ".dat' using 1:3 with lines title label2 lt 2 "<< std::endl;
    plt << "" << std::endl;

    // Data file
    of.open((fname + ".dat").c_str());

    // Header
    of << "#" << std::setw(fw-3) << "X  " << std::setw(fw-3) <<  "Y1" << std::setw(fw-3) <<  "Y2" << std::endl;
    of << std::scientific << std::setprecision(fw-10);
  }

  void add(long double x, long double y1, long double y2)
  {
    of << std::setw(fw) << x << std::setw(fw) << y1 << std::setw(fw) << y2 << std::endl;
  }

}; 

class Plot3
{
  std::string fname;
  std::ofstream  of;
  static const size_t fw = 20;         // Field width

public:
  // fname - w/o extension
  Plot3(const std::string& fname_, const std::string& title, const std::string& xlbl, const std::string& ylbl, 
        const std::string& legend1, const std::string& legend2, const std::string& legend3) : fname(fname_)
  {
    // Plot script
    std::ofstream plt((fname + ".plt").c_str(), std::ofstream::out);
    plt << "#" << std::endl;
    plt << "# Run with gnuplot -p " << fname << ".plt" << std::endl;
    plt << "#" << std::endl;
    plt << "set multiplot" << std::endl;
 
    // plt << "set format y \"%.1f\"" << std::endl;
    plt << "set xtics  norangelimit" << std::endl;
    plt << "set grid x y" << std::endl;
    plt << "set mxtics" << std::endl;
    plt << "set mytics" << std::endl;
    plt << "set title \""<< title << "\"" << std::endl;
    plt << "set xlabel \"" << xlbl << "\"" << std::endl;
    plt << "set ylabel \"" << ylbl << "\"" << std::endl;
    plt << "label1 = \"" << legend1 << "\"" << std::endl;
    plt << "label2 = \"" << legend2 << "\"" << std::endl;
    plt << "label3 = \"" << legend3 << "\"" << std::endl;
    plt << "plot '" << fname << ".dat' using 1:2 with lines title label1 lt 1," 
        << "     '" << fname << ".dat' using 1:3 with lines title label2 lt 2,"
        << "     '" << fname << ".dat' using 1:4 with lines title label3 lt 3 "<< std::endl;
    plt << "" << std::endl;

    // Data file
    of.open((fname + ".dat").c_str());

    // Header
    of << "#" << std::setw(fw-3) << "X  " << std::setw(fw-3) <<  "Y1" << std::setw(fw-3) <<  "Y2" << std::setw(fw-3) <<  "Y3" << std::endl;
    of << std::scientific << std::setprecision(fw-10);
  }

  void add(long double x, long double y1, long double y2, long double y3)
  {
    of << std::setw(fw) << x << std::setw(fw) << y1 << std::setw(fw) << y2 << std::setw(fw) << y3 << std::endl;
  }

}; 

class Plot4
{
  std::string fname;
  std::ofstream  of;
  static const size_t fw = 20;         // Field width

public:
  // fname - w/o extension
  Plot4(const std::string& fname_, const std::string& title, const std::string& xlbl, const std::string& ylbl, 
        const std::string& legend1, const std::string& legend2, const std::string& legend3, const std::string& legend4) : fname(fname_)
  {
    // Plot script
    std::ofstream plt((fname + ".plt").c_str(), std::ofstream::out);
    plt << "#" << std::endl;
    plt << "# Run with gnuplot -p " << fname << ".plt" << std::endl;
    plt << "#" << std::endl;
    plt << "set multiplot" << std::endl;
 
    // plt << "set format y \"%.1f\"" << std::endl;
    plt << "set xtics  norangelimit" << std::endl;
    plt << "set grid x y" << std::endl;
    plt << "set mxtics" << std::endl;
    plt << "set mytics" << std::endl;
    plt << "set title \""<< title << "\"" << std::endl;
    plt << "set xlabel \"" << xlbl << "\"" << std::endl;
    plt << "set ylabel \"" << ylbl << "\"" << std::endl;
    plt << "label1 = \"" << legend1 << "\"" << std::endl;
    plt << "label2 = \"" << legend2 << "\"" << std::endl;
    plt << "label3 = \"" << legend3 << "\"" << std::endl;
    plt << "label4 = \"" << legend4 << "\"" << std::endl;
    plt << "plot '" << fname << ".dat' using 1:2 with lines title label1 lt 1," 
        << "     '" << fname << ".dat' using 1:3 with lines title label2 lt 2,"
        << "     '" << fname << ".dat' using 1:4 with lines title label3 lt 3,"
        << "     '" << fname << ".dat' using 1:5 with lines title label4 lt 4 "<< std::endl;
    plt << "" << std::endl;

    // Data file
    of.open((fname + ".dat").c_str());

    // Header
    of << "#" << std::setw(fw-3) << "X  " << std::setw(fw-3) <<  "Y1" << std::setw(fw-3) <<  "Y2" << std::setw(fw-3) <<  "Y3" << std::setw(fw-3) <<  "Y4" << std::endl;
    of << std::scientific << std::setprecision(fw-10);
  }

  void add(long double x, long double y1, long double y2, long double y3, long double y4)
  {
    of << std::setw(fw) << x << std::setw(fw) << y1 << std::setw(fw) << y2 << std::setw(fw) << y3 << std::setw(fw) << y4 << std::endl;
  }

}; 

class Plot5
{
  std::string fname;
  std::ofstream  of;
  static const size_t fw = 20;         // Field width

public:
  // fname - w/o extension
  Plot5(const std::string& fname_, const std::string& title, const std::string& xlbl, const std::string& ylbl, 
        const std::string& legend1, const std::string& legend2, const std::string& legend3, const std::string& legend4, const std::string& legend5) : fname(fname_)
  {
    // Plot script
    std::ofstream plt((fname + ".plt").c_str(), std::ofstream::out);
    plt << "#" << std::endl;
    plt << "# Run with gnuplot -p " << fname << ".plt" << std::endl;
    plt << "#" << std::endl;
    plt << "set multiplot" << std::endl;
 
    // plt << "set format y \"%.1f\"" << std::endl;
    plt << "set xtics  norangelimit" << std::endl;
    plt << "set grid x y" << std::endl;
    plt << "set mxtics" << std::endl;
    plt << "set mytics" << std::endl;
    plt << "set title \""<< title << "\"" << std::endl;
    plt << "set xlabel \"" << xlbl << "\"" << std::endl;
    plt << "set ylabel \"" << ylbl << "\"" << std::endl;
    plt << "label1 = \"" << legend1 << "\"" << std::endl;
    plt << "label2 = \"" << legend2 << "\"" << std::endl;
    plt << "label3 = \"" << legend3 << "\"" << std::endl;
    plt << "label4 = \"" << legend4 << "\"" << std::endl;
    plt << "label5 = \"" << legend5 << "\"" << std::endl;
    plt << "plot '" << fname << ".dat' using 1:2 with lines title label1 lt 1," 
        << "     '" << fname << ".dat' using 1:3 with lines title label2 lt 2,"
        << "     '" << fname << ".dat' using 1:4 with lines title label3 lt 3,"
        << "     '" << fname << ".dat' using 1:5 with lines title label3 lt 4,"
        << "     '" << fname << ".dat' using 1:6 with lines title label4 lt 5 "<< std::endl;
    plt << "" << std::endl;

    // Data file
    of.open((fname + ".dat").c_str());

    // Header
    of << "#" << std::setw(fw-3) << "X  " << std::setw(fw-3) <<  "Y1" << std::setw(fw-3) <<  "Y2" << std::setw(fw-3) <<  "Y3" << std::setw(fw-3) <<  "Y4" << std::setw(fw-3) <<  "Y5"<< std::endl;
    of << std::scientific << std::setprecision(fw-10);
  }

  void add(long double x, long double y1, long double y2, long double y3, long double y4, long double y5)
  {
    of << std::setw(fw) << x << std::setw(fw) << y1 << std::setw(fw) << y2 << std::setw(fw) << y3 << std::setw(fw) << y4 << std::setw(fw) << y5<< std::endl;
  }

}; 

#endif  // __GNUPLOT_HPP__
