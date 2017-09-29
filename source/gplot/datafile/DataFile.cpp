#include <iostream>
#include <gplot/datafile/DataFile.hpp>

typedef gplot::datafile::DataFile THIS;

THIS::DataFile(std::string filename):
	_M_fp(0),
	_M_filename(filename),
	_M_gate(0)
{
	_M_fp = fopen(filename.c_str(), "w");
}
THIS::~DataFile()
{
	if(_M_fp)
		fclose(_M_fp);
}
void			THIS::pushd(
		unsigned int col,
		double d)
{
	// slot already filled
	assert((_M_gate & (1 << col))==0);

	std::cout << "datafile pushd " << col << " " << d << std::endl;

	_M_gate |= (1 << col);
	
	_M_buf[col] = d;

	write();
}
void			THIS::pushv(
		unsigned int col,
		Eigen::VectorXd d,
		unsigned int col2)
{
	// slot already filled
	assert((_M_gate & (1 << col))==0);

	std::cout << "datafile pushv " << col << " " << d[col2] << " " << col2 << std::endl;

	_M_gate |= (1 << col);
	
	_M_buf[col] = d[col2];

	write();
}
void			THIS::write()
{
	if(_M_gate != ((1u << _M_buf.size()) - 1u)) return;

	_M_gate = 0;

	assert(_M_fp);

	std::cout << "datafile write" << std::endl;

	bool first = true;
	for(unsigned int i = 0; i < _M_buf.size(); ++i)
	{
		if(first)
			first = false;
		else
			fprintf(_M_fp, " ");

		fprintf(_M_fp, "%f", _M_buf[i]);
	}
	fprintf(_M_fp, "\n");

	fflush(_M_fp);

	_M_sig();
}


