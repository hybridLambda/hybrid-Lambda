/* 
 * 
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012, 2013 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-Lambda 
 * 
 * hybrid-Lambda is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
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



