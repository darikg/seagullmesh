from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys, os, glob, re
import setuptools


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

include_dirs = [
    "./include/",
    "./src/docs/",
    # Path to pybind11 headers
    get_pybind_include(),
    get_pybind_include(user=True),
]

library_dirs = None
cgal_libs = ["CGAL", "CGAL_Core"]

conda_prefix = os.getenv('CONDA_PREFIX')
if not conda_prefix:
    conda_prefix = os.getenv('MINICONDAPATH')

if conda_prefix:
    cgal_include = os.path.join(conda_prefix, 'include', 'CGAL')

    if not os.path.exists(cgal_include):
        cgal_include = os.path.join(conda_prefix, 'Library', 'include', 'CGAL')
        include_dirs.append(os.path.join(conda_prefix, 'Library', 'include'))

    eigen_include = os.path.join(conda_prefix, 'include', 'eigen3')
    if not os.path.exists(eigen_include):
        eigen_include = os.path.join(conda_prefix, 'Library', 'include', 'eigen3')

    if os.path.exists(eigen_include):
        print(f"Found eigen in {eigen_include}")
        include_dirs.append(eigen_include)
    else:
        print(f"Eigen not found")

elif os.path.exists('/usr/include/CGAL/'):
    cgal_include = '/usr/include/CGAL/'
else:
    cgal_include = '/usr/local/include/CGAL/'

cgal_version = None
if os.path.exists(os.path.join(cgal_include, 'version.h')):
    with open(os.path.join(cgal_include, 'version.h'), 'r') as f:
        m = re.search(r'#define\s+CGAL_VERSION\s+([\d\.]+)', f.read())
        if m:
            cgal_version = tuple(map(int, m.group(1).split('.')))

if cgal_version:
    print("Found CGAL version: " + '.'.join(map(str, cgal_version)))
else:
    print("Could not determine CGAL version.")
    cgal_version = (5, 0)

if cgal_version >= (5, 0):
    cgal_libs = []  # header only now.

if conda_prefix:

    if sys.platform.startswith('win'):
        prefix = os.path.join(sys.prefix, 'Library\\')
    else:
        prefix = sys.prefix

    print("Looking for CGAL in: ", prefix)

    extra_link_args = []
    if sys.platform == 'darwin' and cgal_version < (5, 0):
        extra_link_args = ['-Wl,-rpath', '-Wl,%s' % os.path.abspath(prefix)]

    if sys.platform == 'win32':
        library_dir = os.path.join(prefix, 'lib')

        if cgal_version < (5, 0):
            # currently we also need to add the apparently windows specific suffix here, it's unclear if this is
            # necessary if CGAL is installed not from CONDA.
            # suffix = "-vc140-mt-4.14.1"
            adjusted_cgal_libs = []
            for lib in cgal_libs:
                candidates = glob.glob(os.path.join(library_dir, lib) + '*.lib')
                if not candidates:
                    raise RuntimeError("Library not found [[{}]]".format(lib))
                c = os.path.basename(candidates[0])
                c = os.path.splitext(c)[0]
                adjusted_cgal_libs.append(c)
            cgal_libs = adjusted_cgal_libs
            print("Names of adjusted CGAL libs: ", adjusted_cgal_libs)

        library_dirs = [library_dir]
        print("Looking for libraries in ", library_dirs)

    if cgal_version < (5, 0):
        include_dirs.insert(1, os.path.join(prefix, 'include'))

if sys.platform == 'darwin' and glob.glob('/usr/local/lib/libboost*-mt*'):
    boost_mt = True
else:
    boost_mt = False

ext_modules = [
    Extension(
        'seagullmesh._seagullmesh',
        [
            'src/util.cpp',
            'src/seagullmesh.cpp',
            'src/mesh.cpp',
            'src/properties.cpp',
            'src/meshing.cpp',
            # 'src/corefine.cpp',
            # 'src/locate.cpp',
            # 'src/parametrize.cpp',
            # 'src/triangulate.cpp',
            # 'src/border.cpp',
            # 'src/simplification.cpp',
            # 'src/skeletonization.cpp',
        ],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        libraries=cgal_libs + ['mpfr',
                   'gmp', 
                   'boost_thread-mt' if boost_mt else 'boost_thread',
                   'boost_atomic-mt' if boost_mt else 'boost_atomic',
                   'boost_system',
                   'boost_date_time',
                   'boost_chrono'],
        language='c++'
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    extra = ['-stdlib=libc++'] if sys.platform == 'darwin' else []
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname] + extra)
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++17',]

    for flag in flags:
        if has_flag(compiler, flag): return flag

    raise RuntimeError('Unsupported compiler -- at least C++17 support '
                       'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc', '/std:c++17'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.14']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def initialize_options(self):
        super(BuildExt, self).initialize_options()
        self.debug = False
        self.parallel = True

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])

        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
            if 'CGAL_DEBUG' not in os.environ:
                opts.append('-DCGAL_DEBUG=1')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
            if 'CGAL_DEBUG' not in os.environ:
                opts.append('/DCGAL_DEBUG=1')

        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts

        build_ext.build_extensions(self)


here = os.path.dirname(os.path.abspath(__file__))
version_ns = {}
with open(os.path.join(here, 'seagullmesh', '_version.py')) as f:
    exec(f.read(), {}, version_ns)

setup(
    name='seagullmesh',
    version=version_ns['__version__'],
    author='Darik Gamble',
    author_email='darik.gamble@gmail.com',
    url='https://github.com/darikg/seagull',
    description="seagullmesh, python bindings to CGAL's surface mesh processing modules",
    long_description='',
    ext_modules=ext_modules,
    install_requires=['pybind11', 'numpy'],
    setup_requires=['pybind11'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    packages=['seagullmesh'],
)
