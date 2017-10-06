#include <cassert>
#include <unistd.h>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <gplot/gplot.hpp>
#include <boost/program_options.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

class Simulation
{
	public:
		virtual void	read(std::string) = 0;
		virtual void	write(std::string) = 0;
		void		run(int frames)
		{

			for(int i = 0; i < frames; ++i)
			{
				calculate(_M_current_frame);
				++_M_current_frame;
			}
		}
		virtual void	calculate(int) = 0;

		int		_M_current_frame;
};

//typedef double FLOAT;
typedef boost::multiprecision::cpp_bin_float_50 FLOAT;

typedef std::complex<FLOAT> C;
typedef Eigen::Matrix<FLOAT,2,1> Vector2F;

class ColorHSV
{
	public:
		ColorHSV(
				FLOAT,
				FLOAT,
				FLOAT);
		FLOAT h, s, v;
};
ColorHSV::ColorHSV(
		FLOAT nh,
		FLOAT ns,
		FLOAT nv):
	h(nh),
	s(ns),
	v(nv)
{
}
class ColorFRGB
{
	public:
		ColorFRGB();
		ColorFRGB(
				FLOAT,
				FLOAT,
				FLOAT);
		ColorFRGB(
				ColorHSV const &);
		ColorFRGB &	operator+=(ColorFRGB const & f);
		FLOAT r, g, b;
};
class ColorCRGB
{
	public:
		ColorCRGB();
		ColorCRGB(unsigned char,unsigned char,unsigned char);
		ColorCRGB(ColorFRGB const &);
		ColorCRGB(ColorHSV const &);
		unsigned char r,g,b;
};
ColorFRGB::ColorFRGB():
	r(0),
	g(0),
	b(0)
{}
ColorFRGB::ColorFRGB(
		FLOAT nr,
		FLOAT ng,
		FLOAT nb):
	r(nr),
	g(ng),
	b(nb)
{
}
ColorCRGB::ColorCRGB():
	r(0),
	g(0),
	b(0)
{
}
ColorFRGB &	ColorFRGB::operator+=(ColorFRGB const & f)
{
	r += f.r;
	g += f.g;
	b += f.b;
	return *this;
}
ColorCRGB::ColorCRGB(
		unsigned char nr, 
		unsigned char ng, 
		unsigned char nb):
	r(nr),
	g(ng),
	b(nb)
{
}
ColorFRGB::ColorFRGB(ColorHSV const & hsv)
{
	r = 0;
	g = 0;
	b = 0;

	FLOAT c = hsv.s * hsv.v;
	FLOAT h = hsv.h * 6.;
	int h1 = (int)h;
	int h2 = h1 % 2;
	FLOAT h3 = h - (h1 - h2);
	
	//FLOAT x1 = (double)std::abs(h2 - 1);
	FLOAT x1 = abs(h3 - (FLOAT)1);

	FLOAT x2 = 1. - x1;
	FLOAT x = c * x2;

	if(h < 1)
	{
		r = c;
		g = x;
	}
	else if(h < 2)
	{
		r = x;
		g = c;
	}
	else if(h < 3)
	{
		g = c;
		b = x;
	}
	else if(h < 4)
	{
		g = x;
		b = c;
	}
	else if(h < 5)
	{
		r = x;
		b = c;
	}
	else
	{
		r = c;
		b = x;
	}
	
	r += hsv.v - c;
	g += hsv.v - c;
	b += hsv.v - c;
}
ColorCRGB::ColorCRGB(ColorHSV const & hsv)
{
	ColorFRGB f(hsv);
	r = (unsigned char)(255 * f.r);
	g = (unsigned char)(255 * f.g);
	b = (unsigned char)(255 * f.b);
}
ColorCRGB::ColorCRGB(ColorFRGB const & f)
{
	r = (unsigned char)(255 * f.r);
	g = (unsigned char)(255 * f.g);
	b = (unsigned char)(255 * f.b);
}


FLOAT			julia_distance(C z, C c)
{
	FLOAT threshold = 1e4;
	FLOAT m;

	C dz(1,0);
	
	for(unsigned int i = 0; i < 1024; ++i)
	{
		dz = z * dz * (FLOAT)2 + (FLOAT)1;
		z = z*z + c;
		m = std::norm(z);

		if(m > threshold) 
		{
			//break;
			return sqrt(m / std::norm(dz)) / 2 * log(m);
		}
	}
	return 0;
}
template<typename T>
T			map(
		T x,
		T x0,
		T x1,
		T y0,
		T y1
		)
{
	T z = (y1 - y0) / (x1 - x0) * (x - x0) + y0;
	z = z > y1 ? y1 : z;
	z = z < y0 ? y0 : z;

	return z;
}

template<typename T> using MatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T> using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;


class CD
{
	public:
		void				resize(int, int);
		void				calculate();
		void				filter();

		void				convert_d_to_c1();
		void				convert_c1_to_c2();

		ColorFRGB			color_func_1(FLOAT d)
		{
			FLOAT logd = log(d);
			FLOAT c = map<FLOAT>(logd, -50, 0, 0, 1);
			return ColorFRGB(c,c,c);
		}
		ColorFRGB			color_func_2(FLOAT d)
		{
			if(d==0) return ColorFRGB();
			FLOAT logd = log(d);
			
			FLOAT x0 = logd_min;
			
			FLOAT hue = map<FLOAT>(logd, x0, 0, 0, 2./3.);
			
			ColorHSV hsv(hue, 1, 1);
			return ColorFRGB(hsv);
		}


		int				_M_sx;
		int				_M_sy;
		Vector2F			_M_center;
		FLOAT				_M_wx;
		
		C				_M_c;
		MatrixX<FLOAT>			_M_d;
		MatrixX<ColorFRGB>		_M_c1;
		MatrixX<ColorCRGB>		_M_c2;

		FLOAT				logd_min;
		FLOAT				logd_max;
		Vector2F			logd_min_coor;

		std::function<ColorFRGB(FLOAT d)>	_M_color_func;

};
void				CD::resize(int sx, int sy)
{
	_M_sx = sx;
	_M_sy = sy;
	_M_d.resize(sx, sy);
	_M_c1.resize(sx, sy);
	_M_c2.resize(sx, sy);
}
ColorFRGB operator/(ColorFRGB const & c, FLOAT f)
{
	return ColorFRGB(
			c.r / f,
			c.g / f,
			c.b / f);
}
void				CD::filter()
{
	int s = 1;
	for(int i = 0; i < _M_sx; ++i)
	{
		for(int j = 0; j < _M_sy; ++j)
		{
			FLOAT d = 0;
			ColorFRGB n;
			n += _M_c1(i,j);
			d += 1;
			for(int i1 = i - s; i1 < (i + s + 1); ++i1)
			{
				if(i1 < 0) continue;
				if(i1 >= _M_sx) continue;

				for(int j1 = j - s; j1 < (j + s + 1); ++j1)
				{
					if(j1 < 0) continue;
					if(j1 >= _M_sy) continue;

					n += _M_c1(i1,j1);
					d += 1;
				}
			}
			_M_c1(i,j) = n / d;
		}
	}

}
void				CD::calculate()
{
	logd_min = std::numeric_limits<FLOAT>::max();
	logd_max = std::numeric_limits<FLOAT>::min();

	FLOAT wy = (FLOAT)_M_sy / (FLOAT)_M_sx * _M_wx;

	FLOAT x, y;
	FLOAT x0 = _M_center[0] - _M_wx/2.;
	FLOAT y0 = _M_center[1] - wy/2.;

	FLOAT dx = _M_wx / (FLOAT)(_M_sx - 1);
	FLOAT dy = wy / (FLOAT)(_M_sy - 1);

	for(int i = 0; i < _M_sx; ++i)
	{
		x = x0 + (FLOAT)i * dx;

		if(0) printf("i = %4u\n", i);

		for(int j = 0; j < _M_sy; ++j)
		{
			y = y0 + (FLOAT)j * dy;

			C z(x, y);
			FLOAT d = julia_distance(z, _M_c);
			FLOAT logd = log(d);

			if(logd < logd_min)
			{
				logd_min = logd;
				logd_min_coor = Vector2F(x,y);
			}
			if(logd > logd_max) logd_max = logd;

			if(0)
			{
				//printf("z = %+6.2f, %+6.2f d = %+8.4f logd = %+8.4f\n", z.real(), z.imag(), d, logd);
			}

			_M_d(i, j) = d;
		}
	}


	std::cout << "wx            " << _M_wx << std::endl;
	std::cout << "logd_min      " << logd_min << std::endl;
	std::cout << "logd_max      " << logd_max << std::endl;
	std::cout << "logd_min_coor " << logd_min_coor << std::endl;

	filter();

	convert_d_to_c1();
	convert_c1_to_c2();
}
void			CD::convert_d_to_c1()
{
	for(int i = 0; i < _M_sx; ++i)
	{
		for(int j = 0; j < _M_sy; ++j)
		{
			FLOAT d = _M_d(i,j);
			//FLOAT logd = log(d);

			assert(_M_color_func);

			_M_c1(i, j) = _M_color_func(d);
		}
	}
}
void			CD::convert_c1_to_c2()
{
	for(int i = 0; i < _M_sx; ++i)
	{
		for(int j = 0; j < _M_sy; ++j)
		{
			_M_c2(i, j) = ColorCRGB(_M_c1(i, j));
		}
	}
}
void			do_plot(
		unsigned int sx,
		unsigned int sy,
		std::shared_ptr<gplot::datafile::DataFile> df,
		MatrixX<ColorCRGB> m,
		unsigned int i)
{
	fseek(df->_M_fp, 0, SEEK_SET);

	df->write(&m(0, 0), sx * sy * sizeof(ColorCRGB));
}
namespace po = boost::program_options;
int main(int ac, char ** av)
{
	std::cout << "float                                     " << sizeof(float) << std::endl;
	std::cout << "double                                    " << sizeof(double) << std::endl;
	std::cout << "boost::multiprecision::cpp_bin_float_50   " << sizeof(boost::multiprecision::cpp_bin_float_50) << std::endl;
	

	CD cd;
	int sx = 100;
	int sy = 100;
	int nframes = 1;

	cd._M_color_func = std::bind(&CD::color_func_2, &cd, std::placeholders::_1);

	//C c(0, 0);
	//C c(-.4, .6);
	//C c(-1, 0);
	cd._M_c = C(-0.70176,-0.3842);

	cd._M_center[0] = 0;//-0.227227;//0.5;
	cd._M_center[1] = 0;//-0.171171;//-0.5;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("f", po::value<int>(&nframes)->default_value(1), "nframes")
		("x", po::value<FLOAT>(&cd._M_center[0]), "x")
		("y", po::value<FLOAT>(&cd._M_center[1]), "y")
		("sx", po::value<int>(&sx)->default_value(100), "sx")
		("sy", po::value<int>(&sy)->default_value(100), "sy")
		("wx", po::value<FLOAT>(&cd._M_wx)->default_value(2), "wx");

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}

	if (vm.count("compression")) {
		std::cout << "Compression level was set to " 
			<< vm["compression"].as<int>() << ".\n";
	} else {
		std::cout << "Compression level was not set.\n";
	}


	cd.resize(sx, sy);

	// datafile
	std::shared_ptr<gplot::datafile::DataFile> df(new gplot::datafile::DataFile("build/test1.bin"));


	// plots	
	std::shared_ptr<gplot::plot::Binary> plot(new gplot::plot::Binary);
	plot->_M_linetype = "rgbimage";
	plot->_M_size_x = sx;
	plot->_M_size_y = sy;

	plot->connect(df, 1, 2);

	// gplot
	gplot::GPlot gp;
	gp.connect(plot);


	std::vector<MatrixX<ColorCRGB>> frames;
	frames.resize(nframes);

	int s = nframes * sx * sy * sizeof(ColorCRGB);
	printf("size: %u mb\n", s/1024/1024);


	for(int i = 0; i < nframes; ++i)
	{
		printf("frame %4u\n", i);

		cd.calculate();

		frames[i] = cd._M_c2;

		cd._M_center = cd.logd_min_coor;
		cd._M_wx *= 0.2;
	}

	while(true)
	{
		for(int i = 0; i < nframes; ++i)
		{
			char buf[100];
			sprintf(buf, "frame %4u/%4u", i, nframes);
			plot->_M_title = buf;
			do_plot(sx, sy, df, frames[i], i);
			getchar();
			//usleep(100000);
		}
	}

	//gp.refresh();

	if(ac > 1) getchar();
}

