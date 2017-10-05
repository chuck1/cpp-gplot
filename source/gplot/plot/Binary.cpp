#include <gplot/datafile/DataFile.hpp>
#include <gplot/GPlot.hpp>
#include <gplot/plot/Binary.hpp>

typedef gplot::plot::Binary THIS;

THIS::Binary()
{
}
void			THIS::plot(
		gplot::GPlot & gp)
{
	assert(_M_file);
	
	gp.write("'%s' binary array=(%u,%u) format='%%uchar' with %s",
			_M_file->_M_filename.c_str(),
			_M_size_x,
			_M_size_y,
			_M_linetype.c_str());

	if(!_M_title.empty())
		gp.write(" title '%s'", _M_title.c_str());
}

