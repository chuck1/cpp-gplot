#include <iostream>
#include <gplot/vector/Vector.hpp>
#include <gplot/datafile/DataFile.hpp>

typedef gplot::datafile::DataFile THIS;

THIS::DataFile(std::string filename):
	_M_fp(0),
	_M_filename(filename)
{
	_M_fp = fopen(filename.c_str(), "w");
}
THIS::~DataFile()
{
	if(_M_fp)
		fclose(_M_fp);
}
void		THIS::write(void* buf, size_t count)
{
	fwrite(buf, count, 1, _M_fp);
	fflush(_M_fp);
	_M_sig();
}

