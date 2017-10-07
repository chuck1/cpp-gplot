#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <ncurses.h>

#include <eigen3/Eigen/Core>
#include <gplot/gplot.hpp>
#include <boost/program_options.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include <color.hpp>

class Interupt
{};

class Simulation
{
	public:
		Simulation():
			_M_current_frame(0) {}

		virtual void	read(std::ifstream &) = 0;
		virtual void	write(std::ofstream &) = 0;

		virtual int	size_frames() = 0;
		virtual void	resize_frames(int) = 0;

		virtual void	calculate(int frame) = 0;

		void		write()
		{
			assert(!_M_filename.empty());
			write(_M_filename);
		}
		void		write(std::string f)
		{
			std::ofstream ofs(f, std::ios::binary);
			write(ofs);
		}
		void		read(std::string f)
		{
			_M_filename = f;
			std::ifstream ifs(f, std::ios::binary);
			read(ifs);
		}
		void		sigint()
		{
			throw Interupt();
		}
		void		run(int frames)
		{
			resize_frames(size_frames() + frames);

			try
			{
				for(int i = 0; i < frames; ++i)
				{
					std::cout << "frame " << _M_current_frame << std::endl;

					calculate(_M_current_frame);

					++_M_current_frame;
				}

			}
			catch(Interupt &)
			{
				std::cout << "interupt caught" << std::endl;
			}

			write();
		}

		int			_M_current_frame;
		std::string		_M_filename;
		static Simulation *	_S_simulation;
};
Simulation * Simulation::_S_simulation;

void intHandler(int)
{
	Simulation::_S_simulation->sigint();
}


//typedef double FLOAT;
typedef boost::multiprecision::cpp_bin_float_50 FLOAT;

typedef std::complex<FLOAT> C;
typedef Eigen::Matrix<FLOAT,2,1> Vector2F;

/*
ColorFRGB operator/(ColorFRGB const & c, FLOAT f)
{
	return ColorFRGB(
			c.r / f,
			c.g / f,
			c.b / f);
}
*/

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


class CD: public Simulation
{
	public:
		using Simulation::read;
		virtual void			read(std::ifstream & ifs)
		{
			ifs.read((char*)&_M_current_frame, sizeof(int));
			ifs.read((char*)&_M_c, sizeof(C));
			ifs.read((char*)&_M_sx, sizeof(int));
			ifs.read((char*)&_M_sy, sizeof(int));
			ifs.read((char*)&_M_center, sizeof(Vector2F));
			ifs.read((char*)&_M_wx, sizeof(FLOAT));
			
			unsigned int s;
			ifs.read((char*)&s, sizeof(unsigned int));
			
			_M_frames.resize(s);
			
			for(unsigned int i = 0; i < s; ++i)
			{
				_M_frames[i].resize(_M_sx, _M_sy);
				ifs.read((char*)&_M_frames[i](0,0), sizeof(ColorCRGB) * _M_sx * _M_sy);
			}

			std::cout << "read" << std::endl;
			info();
		}
		virtual void			write(std::ofstream & ofs)
		{
			ofs.write((char*)&_M_current_frame, sizeof(int));
			ofs.write((char*)&_M_c, sizeof(C));
			ofs.write((char*)&_M_sx, sizeof(int));
			ofs.write((char*)&_M_sy, sizeof(int));
			ofs.write((char*)&_M_center, sizeof(Vector2F));
			ofs.write((char*)&_M_wx, sizeof(FLOAT));
			
			unsigned int s = _M_frames.size();
			ofs.write((char*)&s, sizeof(unsigned int));
			
			for(unsigned int i = 0; i < s; ++i)
			{
				ofs.write((char*)&_M_frames[i](0,0), sizeof(ColorCRGB) * _M_sx * _M_sy);
			}
		}
		virtual int			size_frames()
		{
			return _M_frames.size();
		}
		virtual void			resize_frames(int s)
		{
			_M_frames.resize(s);
		}

		virtual void			calculate(int frame);

		void				info()
		{
			std::cout << "sx            " << _M_sx << std::endl;
			std::cout << "sy            " << _M_sy << std::endl;
			std::cout << "frames size   " << _M_frames.size() << std::endl;
			std::cout << "current frame " << _M_current_frame << std::endl;
		}

		void				resize(int, int);
		void				calculate();
		void				filter();

		void				convert_d_to_c1();
		void				convert_c1_to_c2();

		ColorFRGB<FLOAT>		color_func_1(FLOAT d)
		{
			FLOAT logd = log(d);
			FLOAT c = map<FLOAT>(logd, -50, 0, 0, 1);
			return ColorFRGB<FLOAT>(c,c,c);
		}
		ColorFRGB<FLOAT>		color_func_2(FLOAT d)
		{
			if(d==0) return ColorFRGB<FLOAT>();
			FLOAT logd = log(d);

			FLOAT x0 = logd_min;

			FLOAT hue = map<FLOAT>(logd, x0, 0, 0, 2./3.);

			ColorHSV<FLOAT> hsv(hue, 1, 1);
			return ColorFRGB<FLOAT>(hsv);
		}


		C					_M_c;
		int					_M_sx;
		int					_M_sy;
		Vector2F				_M_center;
		FLOAT					_M_wx;

		std::vector<MatrixX<ColorCRGB>>		_M_frames;


		MatrixX<FLOAT>				_M_d;
		MatrixX<ColorFRGB<FLOAT>>		_M_c1;
		MatrixX<ColorCRGB>			_M_c2;

		FLOAT						logd_min;
		FLOAT						logd_max;
		Vector2F					logd_min_coor;

		std::function<ColorFRGB<FLOAT>(FLOAT d)>	_M_color_func;

};
void				CD::resize(int sx, int sy)
{
	_M_sx = sx;
	_M_sy = sy;
	_M_d.resize(sx, sy);
	_M_c1.resize(sx, sy);
	_M_c2.resize(sx, sy);
}

void				CD::filter()
{
	int s = 1;
	for(int i = 0; i < _M_sx; ++i)
	{
		for(int j = 0; j < _M_sy; ++j)
		{
			FLOAT d = 0;
			ColorFRGB<FLOAT> n;
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
void				CD::calculate(int framei)
{
	logd_min = std::numeric_limits<FLOAT>::max();
	logd_max = std::numeric_limits<FLOAT>::min();

	FLOAT wy = (FLOAT)_M_sy / (FLOAT)_M_sx * _M_wx;

	FLOAT x, y;
	FLOAT x0 = _M_center[0] - _M_wx/2.;
	FLOAT y0 = _M_center[1] - wy/2.;

	FLOAT dx = _M_wx / (FLOAT)(_M_sx - 1);
	FLOAT dy = wy / (FLOAT)(_M_sy - 1);

	printf("\n");
	for(int i = 0; i < _M_sx; ++i)
	{
		x = x0 + (FLOAT)i * dx;


		for(int j = 0; j < _M_sy; ++j)
		{
			mvprintw(0, 0, "i = %6u j = %6u %3.0f%%", i, j, 100.f*(float)(i*_M_sx+j)/(float)(_M_sx*_M_sy));
			refresh();

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
	//printf("\n");


	std::cout << "center        " << std::endl;
	printf("    %+32.24e\n", (double)_M_center[0]);
	printf("    %+32.24e\n", (double)_M_center[1]);
	std::cout << "wx            " << _M_wx << std::endl;
	std::cout << "logd_min      " << logd_min << std::endl;
	std::cout << "logd_max      " << logd_max << std::endl;
	std::cout << "logd_min_coor " << std::endl;
	printf("    %+32.24e\n", (double)logd_min_coor[0]);
	printf("    %+32.24e\n", (double)logd_min_coor[1]);

	filter();

	convert_d_to_c1();
	convert_c1_to_c2();

	_M_frames[framei] = _M_c2;

	// prepare for next frame
	_M_center = logd_min_coor;
	_M_wx *= 0.2;
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
			_M_c2(i, j) = ColorCRGB::convert(_M_c1(i, j));
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
inline bool exists (const std::string& name)
{
	struct stat buffer;   
	return (stat (name.c_str(), &buffer) == 0); 
}
namespace po = boost::program_options;
int main(int ac, char ** av)
{
	initscr();

	signal(SIGINT, intHandler);

	std::cout << "float                                     " << sizeof(float) << std::endl;
	std::cout << "double                                    " << sizeof(double) << std::endl;
	std::cout << "boost::multiprecision::cpp_bin_float_50   " << sizeof(boost::multiprecision::cpp_bin_float_50) << std::endl;

	CD cd;
	Simulation::_S_simulation = &cd;
	cd._M_filename = "build/cd_test.bin";
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
		("file", po::value<std::string>(), "file")
		("f", po::value<int>(&nframes)->default_value(1), "nframes")
		("x", po::value<FLOAT>(&cd._M_center[0]), "x")
		("y", po::value<FLOAT>(&cd._M_center[1]), "y")
		("sx", po::value<int>(&cd._M_sx)->default_value(100), "sx")
		("sy", po::value<int>(&cd._M_sy)->default_value(100), "sy")
		("wx", po::value<FLOAT>(&cd._M_wx)->default_value(2), "wx");

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}


	if(vm.count("file"))
	{
		std::string f = vm["file"].as<std::string>();
		if(exists(f))
			cd.read(f);
		else
			cd._M_filename = f;
	}

	cd.resize(cd._M_sx, cd._M_sy);


	// datafile
	std::shared_ptr<gplot::datafile::DataFile> df(new gplot::datafile::DataFile("build/test1.bin"));


	// plots	
	std::shared_ptr<gplot::plot::Binary> plot(new gplot::plot::Binary);
	plot->_M_linetype = "rgbimage";
	plot->_M_size_x = cd._M_sx;
	plot->_M_size_y = cd._M_sy;

	plot->connect(df, 1, 2);

	// gplot
	gplot::GPlot gp;
	gp.connect(plot);

	int s = nframes * cd._M_sx * cd._M_sy * sizeof(ColorCRGB);
	printf("size: %u mb\n", s/1024/1024);

	cd.run(nframes);

	std::cout << "plot" << std::endl;
	cd.info();

	while(true)
	{
		for(int i = 0; i < cd._M_current_frame; ++i)
		{
			char buf[100];
			sprintf(buf, "frame %4u/%4u", i, cd._M_current_frame);
			plot->_M_title = buf;
			do_plot(cd._M_sx, cd._M_sy, df, cd._M_frames[i], i);
			getchar();
			//usleep(100000);
		}
	}

	//gp.refresh();

	if(ac > 1) getchar();

	endwin();
}

