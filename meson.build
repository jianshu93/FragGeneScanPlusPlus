project('FragGeneScan++', 'c',
  version: '0.1',
  license: 'GPL3+',
  meson_version: '>= 0.38.0',
)

# Dependencies
libmath = meson.get_compiler('c').find_library('m')
threads = dependency('threads')

doxygen = find_program('doxygen', required: false)
if not doxygen.found()
  message('Doxygen could not be found, won\'t build documentation.')
endif
astyle = find_program('astyle', required: false)

# Configuration
conf = configuration_data()
pkgdatadir = join_paths(get_option('prefix'), get_option('datadir'), meson.project_name())

conf.set_quoted('PACKAGE_NAME', meson.project_name())
conf.set_quoted('PACKAGE_STRING', '@0@ @1@'.format(meson.project_name(), meson.project_version()))
conf.set_quoted('PACKAGE_VERSION', meson.project_version())
conf.set_quoted('PKGDATADIR', pkgdatadir)
conf.set_quoted('TOP_SRCDIR', meson.source_root())

# Create a config.h file
configure_file(
  output: 'config.h',
  configuration: conf,
)

# FGS+ executable
fgs_sources = [
  'fasta.c',
  'hmm_lib.c',
  'run_hmm.c',
  'util_lib.c',
]

fgs_flags = [
  '-O3',
  '-march=native',
  '-include', 'config.h',
]

fgs_dependencies = [
  libmath,
  threads,
]

executable('FGS++', fgs_sources,
  dependencies: fgs_dependencies,
  c_args: fgs_flags,
  include_directories: include_directories('.'), # include config.h
  install: true,
)

# Install the train dir
install_subdir('train',
  install_dir: pkgdatadir,
)

# Documentation
if doxygen.found()
  doxyfile = configure_file(
    input: 'Doxyfile.in',
    output: 'Doxyfile',
    configuration: conf,
  )
  custom_target('docs',
    input: doxyfile,
    output: 'doc',
    build_by_default: doxygen.found(),
    command: doxygen,
  )
endif

# Coding conventions
if astyle.found()
  custom_target('astyle',
    input: fgs_sources,
    output: fgs_sources,
    command: [
      astyle,
      '--mode=c',
      '--style=google',
      '-s4',
      '-n',
      '-H',
      '-k3',
      '@INPUT@',
    ],
  )
endif
