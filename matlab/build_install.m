function build_install(path)
if nargin<1
    path = [userpath ,'/SNPLIB'];
end
os = computer('arch');
objext = '.o';
mkdir(path);
linkpath = '';
switch os
    case 'win64'
        optimflags = 'OPTIMFLAGS="/O2 /arch:AVX /Oy- /DNDEBUG "';
        objext = '.obj';        
        ldflags = ' ../openblas/lib/libopenblas.dll.a';
        copyfile('../openblas/bin/*.dll',path);
    case 'glnxa64'
        optimflags = 'CXXOPTIMFLAGS="-std=c++11 -O3 -msse4 -pipe -flto"';
        ldflags = '-lopenblas';
    case 'maci64'
        optimflags = 'CXXOPTIMFLAGS="-std=c++11 -O3 -msse4 -pipe -flto "';
        linkpath = '-L/opt/openblas/lib';
        ldflags = '-lopenblas';
end
filelist = dir('../src/*.cpp');
for i=1:length(filelist)
    mex(optimflags,'-c',['../src/',filelist(i).name]);
end

mex(optimflags,'calc_af_.cpp',['base', objext]);
mex(optimflags,'calc_ibs_.cpp',['barrier', objext],['ibs', objext]);
mex(optimflags,'calc_grm_.cpp',['barrier', objext],['grm', objext]);
mex(optimflags,'calc_hamming_.cpp',['barrier', objext],['hamming',objext]);
mex(optimflags,'calc_pca_loadings_.cpp',['svd',objext],['qr',objext],['pca',objext],linkpath,ldflags);
mex(optimflags,'project_pca_.cpp',['svd',objext],['qr',objext],['pca',objext],linkpath,ldflags);
mex(optimflags,'calc_spectral_ibs_loadings_.cpp',['svd',objext],['qr',objext],['barrier', objext],['ibs',objext],['spectral_ibs',objext],linkpath,ldflags);
mex(optimflags,'project_spectral_ibs_.cpp',['svd',objext],['qr',objext],['barrier', objext],['ibs',objext],['spectral_ibs',objext],linkpath,ldflags);
movefile(['*.',mexext],path);
copyfile('SNPLIB.m',path);
delete(['*',objext]);
cd(path);
end