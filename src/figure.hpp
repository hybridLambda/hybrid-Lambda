#include"net.hpp"

namespace figure{
	class param{
		public:

			param();
			param(int argc, char *argv[]);
			
			void plot(string net_str);
		
		private:
		
			bool plot_bool;
			bool dot_bool;
			int plot_option;//=0;
			bool plot_label;
			bool plot_branch;
			string tex_fig_name;
			string dot_fig_name;
			void set_plot_option_();
			valarray <int>  det_x_node (Net net_dummy);
			void plot_in_latex_(const char* file_name,Net net_dummy);
			void plot_in_latex_file_(Net net_dummy);
			void plot_in_dot_(Net net_dummy);
	};
}



