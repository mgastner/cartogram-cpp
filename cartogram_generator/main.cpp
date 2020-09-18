#include <boost/program_options.hpp>
#include <iostream>

using namespace std;

void on_geometry(std::string geometry_file_name)
{
  cerr << "Using geometry from file " << geometry_file_name << endl;
  return;
}
void on_visual_variables(std::string geometry_file_name)
{
  cerr << "Using visual variables from file " << geometry_file_name << endl;
  return;
}

int main(int argc, const char *argv[])
{
  using namespace boost::program_options;

  // Parse command line options. See
  // https://theboostcpplibraries.com/boost.program_options
  try {
    options_description desc{"Options"};
    desc.add_options()(
      "help,h", "Help screen"
      )(
      "geometry,g",
      value<string>()->notifier(on_geometry),
      "GeoJSON file"
      )(
      "visual_variables,v",
      value<string>()->notifier(on_visual_variables),
      "CSV file with area and (optionally) colour");
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);  // Triggers notifier functions such as on_geometry().
    if (vm.count("help")) {
      std::cout << desc << '\n';
    }
  } catch (const error &ex) {
    std::cerr << ex.what() << '\n';
  }
  return 0;
}
