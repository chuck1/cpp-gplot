#include <cassert>
#include <unistd.h>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <gplot/gplot.hpp>

typedef double FLOAT;

struct Color
{
	unsigned char r,g,b;
};

typedef std::complex<FLOAT> C;

unsigned int		julia(C z, C c)
{
	FLOAT threshold = 100;

	for(unsigned int i = 0; i < 255; ++i)
	{
		if(std::abs(z) > threshold) 
		{
			return i;
		}
		z = z*z + c;
	}
	return 255;
}

Color			f1(FLOAT x, FLOAT y, C c)
{
	unsigned char r = 255 * (sin(x) * 0.5 + 0.5);
	unsigned char g = 255 * (sin(y) * 0.5 + 0.5);
	unsigned char b;

	r = julia(C(x, y), c);
	g = r;
	b = r;

	return Color{r, g, b};
}
Eigen::MatrixXd		make_matrix(
		Eigen::VectorXd x,
		Eigen::VectorXd y,
		std::function<FLOAT(FLOAT, FLOAT)> f
		)
{
	Eigen::MatrixXd m(x.size() + 1, y.size() + 1);
	m(0,0) = y.size();
	for(unsigned int i = 0; i < x.size(); ++i)
		m(i + 1, 0) = x[i];
	for(unsigned int j = 0; j < x.size(); ++j)
		m(0, j + 1) = x[j];
	for(unsigned int i = 0; i < x.size(); ++i)
		for(unsigned int j = 0; j < x.size(); ++j)
			m(i + 1, j + 1) = f(x[i], y[j]);
	return m;
}

template<typename T> using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T> using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename X, typename Y, typename Z>
MatrixX<Z>		make_matrix_uniform(
		VectorX<X> x,
		VectorX<Y> y,
		std::function<Z(X, Y)> f
		)
{
	MatrixX<Z> m(x.size(), y.size());
	for(unsigned int i = 0; i < x.size(); ++i)
	{
		//printf("i = %4u\n", i);
		for(unsigned int j = 0; j < y.size(); ++j)
			m(i, j) = f(x[i], y[j]);
	}
	return m;
}

MatrixX<Color>		do_calc(
		unsigned int sx,
		unsigned int sy,
		FLOAT x_center, 
		FLOAT y_center,
		FLOAT wx)
{
	FLOAT wy = (FLOAT)sy / (FLOAT)sx * wx;

	auto x = VectorX<FLOAT>::LinSpaced(sx, x_center - wx/2, x_center + wx/2);
	auto y = VectorX<FLOAT>::LinSpaced(sy, y_center - wy/2, y_center + wy/2);

	auto m = make_matrix_uniform<FLOAT, FLOAT, Color>(x, y, std::bind(&f1, std::placeholders::_1, std::placeholders::_2, C(-.4, .6)));
	return m;
}
void			do_plot(
		unsigned int sx,
		unsigned int sy,
		std::shared_ptr<gplot::datafile::DataFile> df,
		MatrixX<Color> m,
		unsigned int i)
{
	fseek(df->_M_fp, 0, SEEK_SET);

	df->write(&m(0, 0), sx * sy * sizeof(Color));
}
int main(int ac, char ** av)
{
	//gplot::vector::Vector<FLOAT> v1;
	//gplot::vector::Vector<Eigen::VectorXd> v2;
	
	// datafile
	std::shared_ptr<gplot::datafile::DataFile> df(new gplot::datafile::DataFile("build/test1.bin"));
	//std::shared_ptr<gplot::datafile::DataFile> df2(new gplot::datafile::DataFile("build/test2.bin"));
	
	//df1->connectd(v1, 0);
	//df1->connectv(v2, 1, 0);
	//df1->connectv(v2, 2, 1);
	
	//df2->connectv2(v2, 2);

	unsigned int sx = atoi(av[1]);
	unsigned int sy = atoi(av[2]);
	unsigned int nframes = atoi(av[3]);

	// plots	
	std::shared_ptr<gplot::plot::Binary> plot(new gplot::plot::Binary);
	plot->_M_linetype = "rgbimage";
	plot->_M_size_x = sx;
	plot->_M_size_y = sy;

	//std::shared_ptr<gplot::plot::Plot> p2(new gplot::plot::Plot);
	//std::shared_ptr<gplot::plot::Plot> p3(new gplot::plot::Plot);
	
	//df._M_sig.connect(std::bind(&gplot::plot::Plot::refresh, p1.get()));
	//df._M_sig.connect(std::bind(&gplot::plot::Plot::refresh, p2.get()));
	
	plot->connect(df, 1, 2);
	//p2->connect(df1, 1, 3);
	//p3->connect(df2, 1, 2);
	
	// gplot
	gplot::GPlot gp;

	gp.connect(plot);
	
	//Eigen::MatrixXd m(5,5);
	
	FLOAT x = 0.5;
	FLOAT y = -0.5;
	FLOAT wx = 2;


	std::vector<MatrixX<Color>> frames;
	frames.resize(nframes);

	unsigned int s = nframes * sx * sy * sizeof(Color);
	printf("size: %u mb\n", s/1024/1024);

	for(unsigned int i = 0; i < nframes; ++i)
	{
		printf("frame %4u\n", i);
		frames[i] = do_calc(sx, sy, x, y, wx);

		wx *= 31./32.;
	}

	while(true)
	{
		for(unsigned int i = 0; i < nframes; ++i)
		{
			char buf[100];
			sprintf(buf, "frame %4u/%4u", i, nframes);
			plot->_M_title = buf;
			do_plot(sx, sy, df, frames[i], i);
			usleep(100000);
		}
	}

	//gp.refresh();

	if(ac > 1) getchar();
}

