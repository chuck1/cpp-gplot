#include <iostream>
#include <gplot/plot/Plot.hpp>
#include <gplot/GPlot.hpp>

typedef gplot::GPlot THIS;

THIS::GPlot():
	_M_flags(0),
	_M_fp(0)
{
}
THIS::~GPlot()
{
	if(_M_fp)
		pclose(_M_fp);
}
void			THIS::refresh()
{
	std::cout << "gplot refresh" << std::endl;
	
	if(!_M_fp)
		_M_fp = popen("gnuplot", "w");

	bool first = true;
	for(auto p : _M_plots)
	{
		if(first)
		{
			first = false;
			fprintf(_M_fp, "plot ");
		}
		else
			fprintf(_M_fp, ",");

		p->plot(*this);
	}
	fprintf(_M_fp, "\n");
	fflush(_M_fp);
}


