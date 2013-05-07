//seg-site.hpp

#include"utility.hpp"

namespace seg{
	class param{
		public:
			param();
			param(int argc, char *argv[]);			
			void create_site_data_dir(vector <string> mt_tree_str_s);
			string seg_dir_name;
		private:
			void create_new_site_data(string gt_string_mut_num,int site_i);
			
	};	
	
}
