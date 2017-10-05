#include <iostream>
#include <gplot/vector/Vector.hpp>
#include <gplot/datafile/DataFile2D.hpp>

typedef gplot::datafile::DataFile2D THIS;

THIS::DataFile2D(std::string filename):
	DataFile(filename),
	_M_gate(0)
{
}
THIS::~DataFile2D()
{
}
void			THIS::pushd(
		unsigned int col,
		double d)
{
	// slot already filled
	assert((_M_gate & (1 << col))==0);

	//std::cout << "datafile pushd " << col << " " << d << std::endl;

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

	//std::cout << "datafile pushv " << col << " " << d[col2] << " " << col2 << std::endl;

	_M_gate |= (1 << col);
	
	_M_buf[col] = d[col2];

	write();
}
void			THIS::pushv2(
		Eigen::VectorXd d
		)
{
	// slot already filled
	assert(_M_gate == 0);

	//std::cout << "datafile pushv2 " << std::endl;

	_M_gate = (1 << d.size()) - 1;
	
	_M_buf = d;

	write();
}
void			THIS::write()
{
	if(_M_gate != ((1u << _M_buf.size()) - 1u)) return;

	_M_gate = 0;

	assert(_M_fp);

	//std::cout << "datafile write" << std::endl;

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
void			THIS::connectd(
		gplot::vector::Vector<double> & v,
		unsigned int col
		)
{
	if((col + 1) > _M_buf.size())
		_M_buf.resize(col + 1);

	v._M_sig.connect(std::bind(&gplot::datafile::DataFile2D::pushd, this, col, std::placeholders::_1));
}
void			THIS::connectv(
		gplot::vector::Vector<Eigen::VectorXd> & v,
		unsigned int col,
		unsigned int col1
		)
{
	if((col + 1) > _M_buf.size())
		_M_buf.resize(col + 1);

	v._M_sig.connect(std::bind(&gplot::datafile::DataFile2D::pushv, this, col, std::placeholders::_1, col1));
}
void			THIS::connectv2(
		gplot::vector::Vector<Eigen::VectorXd> & v,
		unsigned int sz
		)
{
	_M_buf.resize(sz);

	v._M_sig.connect(std::bind(&gplot::datafile::DataFile2D::pushv2, this, std::placeholders::_1));
}


