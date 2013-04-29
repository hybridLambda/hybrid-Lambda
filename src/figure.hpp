//figure.hpp

#include"net.hpp"

namespace figure{
	class param{
		public:

		param();
		param(int argc, char *argv[]);
		
		void append_to_log(string logName);

		void plot(string net_str);
		
		private:
		
		bool plot_bool;
		bool dot_bool;
		int plot_option;//=0;
		bool plot_label;
		bool plot_branch;
		string tex_fig_name;
		string dot_fig_name;
		
		//int set_plot_option_(bool plot_label,bool plot_branch);
		void set_plot_option_();

		valarray <int>  det_x_node (Net net_dummy);
		
		void plot_in_latex(const char* file_name, Net net_dummy, int plot_option);
		void plot_in_latex_file(const char* file_name, Net net_dummy, int plot_option);
		void plot_in_dot(const char* file_name, Net net_dummy, int plot_option);


	};
}



