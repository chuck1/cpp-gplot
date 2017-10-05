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
	std::shared_ptr<gplot::datafile::DataFile2D> df1(new gplot::datafile::DataFile2D("build/test1.bin"));
	std::shared_ptr<gplot::datafile::DataFile2D> df2(new gplot::datafile::DataFile2D("build/test2.bin"));
	
	df1->connectd(v1, 0);
	df1->connectv(v2, 1, 0);
	df1->connectv(v2, 2, 1);
	
	df2->connectv2(v2, 2);

	// plots	
	std::shared_ptr<gplot::plot::Plot> p1(new gplot::plot::Plot);
	std::shared_ptr<gplot::plot::Plot> p2(new gplot::plot::Plot);
	std::shared_ptr<gplot::plot::Plot> p3(new gplot::plot::Plot);
	
	//df._M_sig.connect(std::bind(&gplot::plot::Plot::refresh, p1.get()));
	//df._M_sig.connect(std::bind(&gplot::plot::Plot::refresh, p2.get()));
	
	p1->connect(df1, 1, 2);
	p2->connect(df1, 1, 3);
	p3->connect(df2, 1, 2);
	
	// gplot
	gplot::GPlot gp1;
	gplot::GPlot gp2;

	gp1.connect(p1);
	gp1.connect(p2);
	gp2.connect(p3);

	// test
	for(unsigned int i = 0; i < 100; ++i)
	{
		double x = 0 + 0.01 * i;
		v1.push_back(x);
		Eigen::VectorXd v(2); v << sin(x), cos(x);
		v2.push_back(v);
		//usleep(500000);
	}
	
	if(ac > 1) getchar();
}

