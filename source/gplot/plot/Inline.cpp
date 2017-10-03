#include <gplot/datafile/DataFile.hpp>
#include <gplot/GPlot.hpp>
#include <gplot/plot/Inline.hpp>

typedef gplot::plot::Inline THIS;

THIS::Inline()
{
	_M_col_x = 1;
	_M_col_y = 2;
}
void			THIS::plot(
		gplot::GPlot & gp)
{
	gp.write("'-' using %i:%i with %s",
			_M_col_x,
			_M_col_y,
			_M_linetype.c_str());

	if(!_M_title.empty())
		gp.write(" title '%s'", _M_title.c_str());
}
void			THIS::write_data(
		gplot::GPlot & gp)
{
	for(auto p : _M_data)
	{
		gp.write("%f %f\n", p[0], p[1]);
	}
	gp.write("e\n");
}

