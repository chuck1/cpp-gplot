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
void			THIS::open()
{
	_M_fp = popen("gnuplot", "w");
}
void			THIS::flush()
{
	fflush(_M_fp);
}
void			THIS::refresh()
{
	//std::cout << "gplot refresh" << std::endl;
	
	bool first = true;
	for(auto p : _M_plots)
	{
		if(first)
		{
			first = false;
			write("plot ");
		}
		else
			write(",");

		p->plot(*this);
	}
	fprintf(_M_fp, "\n");

	// for inline data
	for(auto p : _M_plots)
	{
		p->write_data(*this);
	}

	fflush(_M_fp);
}
void			THIS::connect(
		std::shared_ptr<gplot::plot::Plot> p
		)
{
	_M_plots.push_back(p);

	p->_M_sig.connect(std::bind(&gplot::GPlot::refresh, this));
}


