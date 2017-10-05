#include <gplot/datafile/DataFile.hpp>
#include <gplot/GPlot.hpp>
#include <gplot/plot/Plot.hpp>

typedef gplot::plot::Plot THIS;

/*
void			THIS::set_x(std::shared_ptr<gplot::dstream::DataStream> ds)
{
	_M_x = ds;
	ds->_M_sig.connect(SIG_REFRESH::slot_type(&THIS::refresh, this).track_foreign(shared_from_this()));
}
void			THIS::set_y(std::shared_ptr<gplot::dstream::DataStream> ds)
{
	_M_y = ds;
}*/
THIS::Plot():
	_M_linetype("lines")
{
}
void			THIS::plot(gplot::GPlot & gp)
{
	assert(_M_file);
	
	gp.write("'%s' using %i:%i with %s",
			_M_file->_M_filename.c_str(),
			_M_col_x,
			_M_col_y,
			_M_linetype.c_str());

	if(!_M_title.empty())
		gp.write(" title '%s'", _M_title.c_str());
}
void			THIS::write_data(
		gplot::GPlot & gp)
{
}
void			THIS::refresh()
{
	_M_sig();
}
void			THIS::connect(
		std::shared_ptr<gplot::datafile::DataFile> df,
		unsigned int x,
		unsigned int y
		)
{
	_M_file = df;
	_M_col_x = x;
	_M_col_y = y;
	df->_M_sig.connect(std::bind(&gplot::plot::Plot::refresh, shared_from_this()));
}


