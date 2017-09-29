#include <cassert>
#include <unistd.h>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <gplot/gplot.hpp>

int main(int ac, char ** av)
{
	gplot::vector::Vector<double> v1;
	gplot::vector::Vector<Eigen::VectorXd> v2;
	
	// datafile
	auto df = std::make_shared<gplot::datafile::DataFile>("test.bin");
	
	df->_M_buf.resize(3);

	v1._M_sig.connect(std::bind(&gplot::datafile::DataFile::pushd, df.get(), 0, std::placeholders::_1));
	v2._M_sig.connect(std::bind(&gplot::datafile::DataFile::pushv, df.get(), 1, std::placeholders::_1, 0));
	v2._M_sig.connect(std::bind(&gplot::datafile::DataFile::pushv, df.get(), 2, std::placeholders::_1, 1));
	
	// plots	
	std::shared_ptr<gplot::plot::Plot> p1(new gplot::plot::Plot);
	std::shared_ptr<gplot::plot::Plot> p2(new gplot::plot::Plot);
	
	//df._M_sig.connect(std::bind(&gplot::plot::Plot::refresh, p1.get()));
	//df._M_sig.connect(std::bind(&gplot::plot::Plot::refresh, p2.get()));
	
	p1->connect(df, 1, 2);
	p2->connect(df, 1, 3);
	
	// gplot
	gplot::GPlot gp;

	gp._M_plots.push_back(p1);
	gp._M_plots.push_back(p2);
	
	p1->_M_sig.connect(std::bind(&gplot::GPlot::refresh, &gp));
	p2->_M_sig.connect(std::bind(&gplot::GPlot::refresh, &gp));

	// test
	for(unsigned int i = 0; i < 1000; ++i)
	{
		double x = 0 + 0.01 * i;
		v1.push_back(x);
		Eigen::VectorXd v(2); v << sin(x), cos(x);
		v2.push_back(v);
		//usleep(500000);
	}

	getchar();
}

