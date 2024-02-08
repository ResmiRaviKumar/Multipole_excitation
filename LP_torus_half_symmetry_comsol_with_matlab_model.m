% plane wave - torus interaction
% half symmetry model
%% PREAMBLE

%clc; % clear command window
clear all; close all;
tStart=tic;

% Parameters defined for Matlab
q1=1;
small_r=(1/4)*q1;%q1;
R_1=(1/2)*q1;
epsilon0=8.85418782e-12; % m-3 kg-1 s4 A2
n=4;
epsilonr=n.^2; % dielectric constant of the material

%% wavelength
min_lamda=1*q1;%0.5
max_lamda=10*q1;
Nlamda=200;%400;%360;%10 points per wavelength for the shortest wavelength
c=2.998e8;
min_nu=c/max_lamda;
max_nu=c/min_lamda;
dnu=(max_nu-min_nu)/(Nlamda-1);
nnu=min_nu:dnu:max_nu;
llamda=c./nnu;

llamda=llamda(181:200);
% %  llamda=llamda(5);
% llamda=llamda(49:50);
% llamda=1.2727*q1;%2.6432*q1;%1.4157*q1;%10*q1%16.8*q1;%
% c=2.998e8;

% Geometry is adjusted for focused beam
for jlamda=1:numel(llamda)
    lamda=llamda(jlamda);
    lda0=lamda;
    
tic

       r_sph=(small_r+R_1)+3*lda0;
       r_sph_outer=r_sph+0.2*r_sph;
       
       
%          r_sph=(small_r+R_1)+3*lda0;
%        r_sph_outer=r_sph+0.2*r_sph;
                  
% mesh size is decided in such a way that we need atleast 10 points per wavelength for the shortest wavelength and atleast 10 in the smallest small_r
   hmax_air=lamda/10;
   hmin_air=hmax_air/2;
   
    if lamda/10/n<small_r/10
       hmax_torus=lamda/10/n; 
       hmin_torus=hmax_torus/2;
     else
       hmax_torus=small_r/10;
       hmin_torus=hmax_torus/2;
    end 

clear model;

import com.comsol.model.*
import com.comsol.model.util.*
model=ModelUtil.create('Model');

%% Parameters
% all the parameters are normalized to q1
%% Parameters defined for comsol

model.param.set('q1', num2str(q1), 'parameter of wavlength');
model.param.set('R_1',  num2str(R_1), 'torus major radius');
model.param.set('small_r', num2str(small_r), 'torus minor radius');

model.param.set('n', num2str(n), 'r.i of dielectric');
model.param.set('epsilon0', num2str(epsilon0), 'm-3 kg-1 s4 A2');
model.param.set('epsilonr', num2str(epsilonr), 'permittivity of dielectric');
model.param.set('na', '1');

model.param.set('hmin_torus', num2str(hmin_torus));
model.param.set('hmax_torus', num2str(hmax_torus));
model.param.set('hmin_air', num2str(hmin_air));
model.param.set('hmax_air', num2str(hmax_air));

model.param.set('r_sph', num2str(r_sph));
model.param.set('r_sph_outer', num2str(r_sph_outer));

model.param.set('lda0', num2str(lamda));


model.param.set('I0', '1[MW/m^2]', 'Intensity of incident field');
model.param.set('delta1', '1');
model.param.set('delta0', '0');
model.param.set('w', '2*pi*(c_const/lda0)', 'angular frequency');
model.param.set('k', 'w/c_const', 'wave number');

%% Geometry

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);

model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('pos', {'-(r_sph_outer)' '0' '-(r_sph_outer)'});
model.component('comp1').geom('geom1').feature('blk1').set('size', {'2*(r_sph_outer)' '2*(r_sph_outer)' '2*(r_sph_outer)'});
model.component('comp1').geom('geom1').feature('blk1').set('layerbottom', false);

model.component('comp1').geom('geom1').create('sph2', 'Sphere');
model.component('comp1').geom('geom1').feature('sph2').set('layername', {'Layer 1' 'Layer 2'});
model.component('comp1').geom('geom1').feature('sph2').set('layer', {'(0.15*r_sph)' '(0.2*r_sph)'});
model.component('comp1').geom('geom1').feature('sph2').set('r', 'r_sph_outer');

model.component('comp1').geom('geom1').create('int1', 'Intersection');
model.component('comp1').geom('geom1').feature('int1').selection('input').set({'blk1' 'sph2'});

model.component('comp1').geom('geom1').create('blk3', 'Block');
model.component('comp1').geom('geom1').feature('blk3').set('pos', {'-r_sph' '0' '-r_sph'});
model.component('comp1').geom('geom1').feature('blk3').set('size', {'2*r_sph' '2*r_sph' '2*r_sph'});
model.component('comp1').geom('geom1').feature('blk3').set('layerbottom', false);

model.component('comp1').geom('geom1').create('tor1', 'Torus');
model.component('comp1').geom('geom1').feature('tor1').set('rmaj', 'q1/2');
model.component('comp1').geom('geom1').feature('tor1').set('rmin', 'q1/4');

model.component('comp1').geom('geom1').create('int2', 'Intersection');
model.component('comp1').geom('geom1').feature('int2').selection('input').set({'blk3' 'tor1'});
model.component('comp1').geom('geom1').create('uni1', 'Union');
model.component('comp1').geom('geom1').feature('uni1').selection('input').set({'int1'});
model.component('comp1').geom('geom1').run;

% geomcheck=1;
% if geomcheck==1
%    figure; mphgeom(model)
% end

%% Selection
model.component('comp1').selection.create('sel2', 'Explicit');
model.component('comp1').selection('sel2').label('NAnoparticle');
model.component('comp1').selection('sel2').set([6]);

model.component('comp1').selection.create('sel4', 'Explicit');
model.component('comp1').selection('sel4').label('sphere1');
model.component('comp1').selection('sel4').set([5]);

model.component('comp1').selection.create('sel6', 'Explicit');
model.component('comp1').selection('sel6').label('inner_sphere_surface');
model.component('comp1').selection('sel6').geom('geom1', 2);
model.component('comp1').selection('sel6').set([12 13 25 28]);

model.component('comp1').selection.create('sel7', 'Explicit');
model.component('comp1').selection('sel7').label('nanoparticle surface');
model.component('comp1').selection('sel7').geom('geom1', 2);
model.component('comp1').selection('sel7').set([14 15 16 17 18 33 34 35 36 39]);

model.component('comp1').selection.create('sel8', 'Explicit');
model.component('comp1').selection('sel8').label('PML_outer_surface');
model.component('comp1').selection('sel8').geom('geom1', 2);
model.component('comp1').selection('sel8').set([4 5 21 32]);

model.component('comp1').selection.create('sel9', 'Explicit');
model.component('comp1').selection('sel9').label('sphere 2');
model.component('comp1').selection('sel9').set([3 4 8 9]);

model.component('comp1').selection.create('sel10', 'Explicit');
model.component('comp1').selection('sel10').label('PML domain');
model.component('comp1').selection('sel10').set([1 2 7 10]);


model.component('comp1').selection.create('sel11', 'Explicit');
model.component('comp1').selection('sel11').geom('geom1', 2);
model.component('comp1').selection('sel11').set([9 10 24 31]);
model.component('comp1').selection('sel11').label('outer_sphere_surface');


model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.named('sel10');
model.component('comp1').coordSystem('pml1').set('ScalingType', 'Spherical');

%% Mesh
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.named('sel2');
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').selection.set([1 2 3 4 5 7 8 9 10]);
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftet1').set('optcurved', false);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', 'hmax_torus');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', 'hmin_torus');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);

model.component('comp1').mesh('mesh1').feature('ftet2').set('optcurved', false);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmax', 'hmax_air');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmin', 'hmin_air');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').run;

%% Integration

model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').selection.named('sel2');
model.component('comp1').cpl('intop1').set('opname', 'intop_vol');

model.component('comp1').cpl.create('intop2', 'Integration');
model.component('comp1').cpl('intop2').selection.named('sel6');
model.component('comp1').cpl('intop2').set('opname', 'intop_surf');

model.component('comp1').cpl.create('intop3', 'Integration');
model.component('comp1').cpl('intop3').selection.named('sel7');
model.component('comp1').cpl('intop3').set('opname', 'intop_nanosurf');

%% Variables 

model.component('comp1').variable.create('var2');
model.component('comp1').variable('var2').set('epsilon', '(n)^2');
model.component('comp1').variable('var2').set('const', '-1i*w*epsilon0_const*(epsilon-1)');
model.component('comp1').variable('var2').set('Jx', 'const*ewfd.Ex', 'Current density, x component');
model.component('comp1').variable('var2').set('Jy', 'const*ewfd.Ey', 'Current density, y component');
model.component('comp1').variable('var2').set('Jz', 'const*ewfd.Ez', 'Current density, z component');
model.component('comp1').variable('var2').set('conversion', '((c_const^2)/2)*(mu0_const/(4*pi))');
model.component('comp1').variable('var2').label('current density');

model.component('comp1').variable.create('var3');
model.component('comp1').variable('var3').set('px', '(1i/w)*intop_vol(Jx)');
model.component('comp1').variable('var3').set('py', '(1i/w)*intop_vol(Jy)');
model.component('comp1').variable('var3').set('pz', '(1i/w)*intop_vol(Jz)');
model.component('comp1').variable('var3').set('p', 'sqrt(abs(px)^2+abs(py)^2+abs(pz)^2)');
model.component('comp1').variable('var3').label('electric dipole moment');

model.component('comp1').variable.create('var5');
model.component('comp1').variable('var5').set('rdotj', '(x*Jx)+(y*Jy)+(z*Jz)');
model.component('comp1').variable('var5').set('rsqr', 'abs(x^2)+abs(y^2)+abs(z^2)');
model.component('comp1').variable('var5').set('tx', '(rdotj*x)-(2*rsqr*Jx)');
model.component('comp1').variable('var5').set('Tx', '(1/10)*intop_vol(tx)');
model.component('comp1').variable('var5').set('ty', '(rdotj*y)-(2*rsqr*Jy)');
model.component('comp1').variable('var5').set('Ty', '(1/10)*intop_vol(ty)');
model.component('comp1').variable('var5').set('tz', '(rdotj*z)-(2*rsqr*Jz)');
model.component('comp1').variable('var5').set('Tz', '(1/10)*intop_vol(tz)');
model.component('comp1').variable('var5').set('T', 'sqrt(abs(Tx)^2+abs(Ty)^2+abs(Tz)^2)');
model.component('comp1').variable('var5').label('toroidal dipole moment');

model.component('comp1').variable.create('var11');
model.component('comp1').variable('var11').set('tx1', '((3*rsqr*Jx)-(2*rdotj*x))*rsqr');
model.component('comp1').variable('var11').set('Tx1', '(1/280)*intop_vol(tx1)');
model.component('comp1').variable('var11').set('ty1', '((3*rsqr*Jy)-(2*rdotj*y))*rsqr');
model.component('comp1').variable('var11').set('Ty1', '(1/280)*intop_vol(ty1)');
model.component('comp1').variable('var11').set('tz1', '((3*rsqr*Jz)-(2*rdotj*z))*rsqr');
model.component('comp1').variable('var11').set('Tz1', '(1/280)*intop_vol(tz1)');
model.component('comp1').variable('var11').set('T1', 'sqrt(abs(Tx1)^2+abs(Ty1)^2+abs(Tz1)^2)');
model.component('comp1').variable('var11').label('toroidal dipole msr');

model.component('comp1').variable.create('var4');
model.component('comp1').variable('var4').set('m1', '(y*Jz-z*Jy)');
model.component('comp1').variable('var4').set('m2', '(z*Jx-x*Jz)');
model.component('comp1').variable('var4').set('m3', '(x*Jy-y*Jx)');
model.component('comp1').variable('var4').set('mx', '(1/(2))*intop_vol(m1)');
model.component('comp1').variable('var4').set('my', '(1/(2))*intop_vol(m2)');
model.component('comp1').variable('var4').set('mz', '(1/(2))*intop_vol(m3)');
model.component('comp1').variable('var4').set('m', 'sqrt(abs(mx)^2+abs(my)^2+abs(mz)^2)');
model.component('comp1').variable('var4').label('magnetic dipole moment');

model.component('comp1').variable.create('var18');
model.component('comp1').variable('var18').set('tx_m', 'rsqr*m1');
model.component('comp1').variable('var18').set('Tx_m', '((1i*w)/20)*intop_vol(tx_m)');
model.component('comp1').variable('var18').set('ty_m', 'rsqr*m2');
model.component('comp1').variable('var18').set('Ty_m', '((1i*w)/20)*intop_vol(ty_m)');
model.component('comp1').variable('var18').set('tz_m', 'rsqr*m3');
model.component('comp1').variable('var18').set('Tz_m', '((1i*w)/20)*intop_vol(tz_m)');
model.component('comp1').variable('var18').set('T_m', 'sqrt(abs(Tx_m)^2+abs(Ty_m)^2+abs(Tz_m)^2)');
model.component('comp1').variable('var18').label('magnetic toroidal dipole moment');

model.component('comp1').variable.create('var6');
model.component('comp1').variable('var6').set('qxx', '(x*Jx)+(x*Jx)');
model.component('comp1').variable('var6').set('qxx_e', 'qxx-((2/3)*delta1*rdotj)');
model.component('comp1').variable('var6').set('Qxx_e', '(1i/w)*intop_vol(qxx_e)');
model.component('comp1').variable('var6').set('qxy', '(x*Jy)+(y*Jx)');
model.component('comp1').variable('var6').set('qxy_e', 'qxy-((2/3)*delta0*rdotj)');
model.component('comp1').variable('var6').set('Qxy_e', '(1i/w)*intop_vol(qxy_e)');
model.component('comp1').variable('var6').set('qxz', '(x*Jz)+(z*Jx)');
model.component('comp1').variable('var6').set('qxz_e', 'qxz-((2/3)*delta0*rdotj)');
model.component('comp1').variable('var6').set('Qxz_e', '(1i/w)*intop_vol(qxz_e)');
model.component('comp1').variable('var6').set('qyx', '(y*Jx)+(x*Jy)');
model.component('comp1').variable('var6').set('qyx_e', 'qyx-((2/3)*delta0*rdotj)');
model.component('comp1').variable('var6').set('Qyx_e', '(1i/w)*intop_vol(qyx_e)');
model.component('comp1').variable('var6').set('qyy', '(y*Jy)+(y*Jy)');
model.component('comp1').variable('var6').set('qyy_e', 'qyy-((2/3)*delta1*rdotj)');
model.component('comp1').variable('var6').set('Qyy_e', '(1i/w)*intop_vol(qyy_e)');
model.component('comp1').variable('var6').set('qyz', '(y*Jz)+(z*Jy)');
model.component('comp1').variable('var6').set('qyz_e', 'qyz-((2/3)*delta0*rdotj)');
model.component('comp1').variable('var6').set('Qyz_e', '(1i/w)*intop_vol(qyz_e)');
model.component('comp1').variable('var6').set('qzx', '(z*Jx)+(x*Jz)');
model.component('comp1').variable('var6').set('qzx_e', 'qzx-((2/3)*delta0*rdotj)');
model.component('comp1').variable('var6').set('Qzx_e', '(1i/w)*intop_vol(qzx_e)');
model.component('comp1').variable('var6').set('qzy', '(z*Jy)+(y*Jz)');
model.component('comp1').variable('var6').set('qzy_e', 'qzy-((2/3)*delta0*rdotj)');
model.component('comp1').variable('var6').set('Qzy_e', '(1i/w)*intop_vol(qzy_e)');
model.component('comp1').variable('var6').set('qzz', '(z*Jz)+(z*Jz)');
model.component('comp1').variable('var6').set('qzz_e', 'qzz-((2/3)*delta1*rdotj)');
model.component('comp1').variable('var6').set('Qzz_e', '(1i/w)*intop_vol(qzz_e)');
model.component('comp1').variable('var6').set('QQe', 'sqrt(abs(Qxx_e)^2+abs(Qxy_e)^2+abs(Qxz_e)^2+abs(Qyx_e)^2+abs(Qyy_e)^2+abs(Qyz_e)^2+abs(Qzx_e)^2+abs(Qzy_e)^2+abs(Qzz_e)^2)');
model.component('comp1').variable('var6').label('electric quadrupole moment');

model.component('comp1').variable.create('var12');
model.component('comp1').variable('var12').set('qxx_t', '(4*x*x*rdotj)-(5*rsqr*(x*Jx+x*Jx))+(2*rsqr*rdotj*delta1)');
model.component('comp1').variable('var12').set('Qxx_t', '(1/42)*intop_vol(qxx_t)');
model.component('comp1').variable('var12').set('qxy_t', '(4*x*y*rdotj)-(5*rsqr*(x*Jy+y*Jx))+(2*rsqr*rdotj*delta0)');
model.component('comp1').variable('var12').set('Qxy_t', '(1/42)*intop_vol(qxy_t)');
model.component('comp1').variable('var12').set('qxz_t', '(4*x*z*rdotj)-(5*rsqr*(x*Jz+z*Jx))+(2*rsqr*rdotj*delta0)');
model.component('comp1').variable('var12').set('Qxz_t', '(1/42)*intop_vol(qxz_t)');
model.component('comp1').variable('var12').set('qyx_t', '(4*y*x*rdotj)-(5*rsqr*(y*Jx+x*Jy))+(2*rsqr*rdotj*delta0)');
model.component('comp1').variable('var12').set('Qyx_t', '(1/42)*intop_vol(qyx_t)');
model.component('comp1').variable('var12').set('qyy_t', '(4*y*y*rdotj)-(5*rsqr*(y*Jy+y*Jy))+(2*rsqr*rdotj*delta1)');
model.component('comp1').variable('var12').set('Qyy_t', '(1/42)*intop_vol(qyy_t)');
model.component('comp1').variable('var12').set('qyz_t', '(4*y*z*rdotj)-(5*rsqr*(y*Jz+z*Jy))+(2*rsqr*rdotj*delta0)');
model.component('comp1').variable('var12').set('Qyz_t', '(1/42)*intop_vol(qyz_t)');
model.component('comp1').variable('var12').set('qzx_t', '(4*z*x*rdotj)-(5*rsqr*(z*Jx+x*Jz))+(2*rsqr*rdotj*delta0)');
model.component('comp1').variable('var12').set('Qzx_t', '(1/42)*intop_vol(qzx_t)');
model.component('comp1').variable('var12').set('qzy_t', '(4*z*y*rdotj)-(5*rsqr*(z*Jy+y*Jz))+(2*rsqr*rdotj*delta0)');
model.component('comp1').variable('var12').set('Qzy_t', '(1/42)*intop_vol(qzy_t)');
model.component('comp1').variable('var12').set('qzz_t', '(4*z*z*rdotj)-(5*rsqr*(z*Jz+z*Jz))+(2*rsqr*rdotj*delta1)');
model.component('comp1').variable('var12').set('Qzz_t', '(1/42)*intop_vol(qzz_t)');
model.component('comp1').variable('var12').set('QQt', 'sqrt(abs(Qxx_t)^2+abs(Qxy_t)^2+abs(Qxz_t)^2+abs(Qyx_t)^2+abs(Qyy_t)^2+abs(Qyz_t)^2+abs(Qzx_t)^2+abs(Qzy_t)^2+abs(Qzz_t)^2)');
model.component('comp1').variable('var12').label('toroidal quadrupole moment');

model.component('comp1').variable.create('var7');
model.component('comp1').variable('var7').set('qxx_m', 'm1*x+m1*x');
model.component('comp1').variable('var7').set('Qxx_m', '(1/(3))*intop_vol(qxx_m)');
model.component('comp1').variable('var7').set('qxy_m', 'm1*y+m2*x');
model.component('comp1').variable('var7').set('Qxy_m', '(1/(3))*intop_vol(qxy_m)');
model.component('comp1').variable('var7').set('qxz_m', 'm1*z+m3*x');
model.component('comp1').variable('var7').set('Qxz_m', '(1/(3))*intop_vol(qxz_m)');
model.component('comp1').variable('var7').set('qyx_m', 'm2*x+m1*y');
model.component('comp1').variable('var7').set('Qyx_m', '(1/(3))*intop_vol(qyx_m)');
model.component('comp1').variable('var7').set('qyy_m', 'm2*y+m2*y');
model.component('comp1').variable('var7').set('Qyy_m', '(1/(3))*intop_vol(qyy_m)');
model.component('comp1').variable('var7').set('qyz_m', 'm2*z+m3*y');
model.component('comp1').variable('var7').set('Qyz_m', '(1/(3))*intop_vol(qyz_m)');
model.component('comp1').variable('var7').set('qzx_m', 'm3*x+m1*z');
model.component('comp1').variable('var7').set('Qzx_m', '(1/(3))*intop_vol(qzx_m)');
model.component('comp1').variable('var7').set('qzy_m', 'm3*y+m2*z');
model.component('comp1').variable('var7').set('Qzy_m', '(1/(3))*intop_vol(qzy_m)');
model.component('comp1').variable('var7').set('qzz_m', 'm3*z+m3*z');
model.component('comp1').variable('var7').set('Qzz_m', '(1/(3))*intop_vol(qzz_m)');
model.component('comp1').variable('var7').set('QQm', 'sqrt(abs(Qxx_m)^2+abs(Qxy_m)^2+abs(Qxz_m)^2+abs(Qyx_m)^2+abs(Qyy_m)^2+abs(Qyz_m)^2+abs(Qzx_m)^2+abs(Qzy_m)^2+abs(Qzz_m)^2)');
model.component('comp1').variable('var7').label('magnetic quadrupole moment');

model.component('comp1').variable.create('var17');
model.component('comp1').variable('var17').set('qxx_tm', 'rsqr*(m1*x+m1*x)');
model.component('comp1').variable('var17').set('Qxx_tm', '((1i*w)/42)*intop_vol(qxx_tm)');
model.component('comp1').variable('var17').set('qxy_tm', 'rsqr*(m1*y+m2*x)');
model.component('comp1').variable('var17').set('Qxy_tm', '((1i*w)/42)*intop_vol(qxy_tm)');
model.component('comp1').variable('var17').set('qxz_tm', 'rsqr*(m1*z+m3*x)');
model.component('comp1').variable('var17').set('Qxz_tm', '((1i*w)/42)*intop_vol(qxz_tm)');
model.component('comp1').variable('var17').set('qyx_tm', 'rsqr*(m2*x+m1*y)');
model.component('comp1').variable('var17').set('Qyx_tm', '((1i*w)/42)*intop_vol(qyx_tm)');
model.component('comp1').variable('var17').set('qyy_tm', 'rsqr*(m2*y+m2*y)');
model.component('comp1').variable('var17').set('Qyy_tm', '((1i*w)/42)*intop_vol(qyy_tm)');
model.component('comp1').variable('var17').set('qyz_tm', 'rsqr*(m2*z+m3*y)');
model.component('comp1').variable('var17').set('Qyz_tm', '((1i*w)/42)*intop_vol(qyz_tm)');
model.component('comp1').variable('var17').set('qzx_tm', 'rsqr*(m3*x+m1*z)');
model.component('comp1').variable('var17').set('Qzx_tm', '((1i*w)/42)*intop_vol(qzx_tm)');
model.component('comp1').variable('var17').set('qzy_tm', 'rsqr*(m3*y+m2*z)');
model.component('comp1').variable('var17').set('Qzy_tm', '((1i*w)/42)*intop_vol(qzy_tm)');
model.component('comp1').variable('var17').set('qzz_tm', 'rsqr*(m3*z+m3*z)');
model.component('comp1').variable('var17').set('Qzz_tm', '((1i*w)/42)*intop_vol(qzz_tm)');
model.component('comp1').variable('var17').set('QQ_tm', 'sqrt(abs(Qxx_tm)^2+abs(Qxy_tm)^2+abs(Qxz_tm)^2+abs(Qyx_tm)^2+abs(Qyy_tm)^2+abs(Qyz_tm)^2+abs(Qzx_tm)^2+abs(Qzy_tm)^2+abs(Qzz_tm)^2)');
model.component('comp1').variable('var17').label('magnetic toroidal quadrupole');

model.component('comp1').variable.create('var10');
model.component('comp1').variable('var10').set('oxxx_e', 'Jx*(((x*x)/3)-((1/5)*rsqr*delta1))+x*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))+Jx*(((x*x)/3)-((1/5)*rsqr*delta1))+x*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))+Jx*(((x*x)/3)-((1/5)*rsqr*delta1))+x*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Oxxx_e', '(i/w)*intop_vol(oxxx_e)');
model.component('comp1').variable('var10').set('oxxy_e', 'Jx*(((x*y)/3)-((1/5)*rsqr*delta0))+x*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))+Jx*(((x*y)/3)-((1/5)*rsqr*delta0))+x*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))+Jy*(((x*x)/3)-((1/5)*rsqr*delta1))+y*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Oxxy_e', '(i/w)*intop_vol(oxxy_e)');
model.component('comp1').variable('var10').set('oxxz_e', 'Jx*(((x*z)/3)-((1/5)*rsqr*delta0))+x*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))+Jx*(((x*z)/3)-((1/5)*rsqr*delta0))+x*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))+Jz*(((x*x)/3)-((1/5)*rsqr*delta1))+z*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Oxxz_e', '(i/w)*intop_vol(oxxz_e)');
model.component('comp1').variable('var10').set('oxyx_e', 'Jx*(((y*x)/3)-((1/5)*rsqr*delta0))+x*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))+Jy*(((x*x)/3)-((1/5)*rsqr*delta1))+y*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))+Jx*(((y*x)/3)-((1/5)*rsqr*delta0))+x*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oxyx_e', '(i/w)*intop_vol(oxyx_e)');
model.component('comp1').variable('var10').set('oxyy_e', 'Jx*(((y*y)/3)-((1/5)*rsqr*delta1))+x*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))+Jy*(((x*y)/3)-((1/5)*rsqr*delta0))+y*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))+Jy*(((y*x)/3)-((1/5)*rsqr*delta0))+y*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oxyy_e', '(i/w)*intop_vol(oxyy_e)');
model.component('comp1').variable('var10').set('oxyz_e', 'Jx*(((y*z)/3)-((1/5)*rsqr*delta0))+x*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))+Jy*(((x*z)/3)-((1/5)*rsqr*delta0))+y*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))+Jz*(((y*x)/3)-((1/5)*rsqr*delta0))+z*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oxyz_e', '(i/w)*intop_vol(oxyz_e)');
model.component('comp1').variable('var10').set('oxzx_e', 'Jx*(((z*x)/3)-((1/5)*rsqr*delta0))+x*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))+Jz*(((x*x)/3)-((1/5)*rsqr*delta1))+z*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))+Jx*(((z*x)/3)-((1/5)*rsqr*delta0))+x*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oxzx_e', '(i/w)*intop_vol(oxzx_e)');
model.component('comp1').variable('var10').set('oxzy_e', 'Jx*(((z*y)/3)-((1/5)*rsqr*delta0))+x*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))+Jz*(((x*y)/3)-((1/5)*rsqr*delta0))+z*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))+Jy*(((z*x)/3)-((1/5)*rsqr*delta0))+y*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oxzy_e', '(i/w)*intop_vol(oxzy_e)');
model.component('comp1').variable('var10').set('oxzz_e', 'Jx*(((z*z)/3)-((1/5)*rsqr*delta1))+x*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))+Jz*(((x*z)/3)-((1/5)*rsqr*delta0))+z*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))+Jz*(((z*x)/3)-((1/5)*rsqr*delta0))+z*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oxzz_e', '(i/w)*intop_vol(oxzz_e)');
model.component('comp1').variable('var10').set('oyxx_e', 'Jy*(((x*x)/3)-((1/5)*rsqr*delta1))+y*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))+Jx*(((y*x)/3)-((1/5)*rsqr*delta0))+x*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))+Jx*(((x*y)/3)-((1/5)*rsqr*delta0))+x*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oyxx_e', '(i/w)*intop_vol(oyxx_e)');
model.component('comp1').variable('var10').set('oyxy_e', 'Jy*(((x*y)/3)-((1/5)*rsqr*delta0))+y*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))+Jx*(((y*y)/3)-((1/5)*rsqr*delta1))+x*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))+Jy*(((x*y)/3)-((1/5)*rsqr*delta1))+y*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oyxy_e', '(i/w)*intop_vol(oyxy_e)');
model.component('comp1').variable('var10').set('oyxz_e', 'Jy*(((x*z)/3)-((1/5)*rsqr*delta0))+y*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))+Jx*(((y*z)/3)-((1/5)*rsqr*delta0))+x*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))+Jz*(((x*y)/3)-((1/5)*rsqr*delta0))+z*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oyxz_e', '(i/w)*intop_vol(oyxz_e)');
model.component('comp1').variable('var10').set('oyyx_e', 'Jy*(((y*x)/3)-((1/5)*rsqr*delta0))+y*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))+Jy*(((y*x)/3)-((1/5)*rsqr*delta0))+y*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))+Jx*(((y*y)/3)-((1/5)*rsqr*delta1))+x*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Oyyx_e', '(i/w)*intop_vol(oyyx_e)');
model.component('comp1').variable('var10').set('oyyy_e', 'Jy*(((y*y)/3)-((1/5)*rsqr*delta1))+y*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))+Jy*(((y*y)/3)-((1/5)*rsqr*delta1))+y*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))+Jy*(((y*y)/3)-((1/5)*rsqr*delta1))+y*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Oyyy_e', '(i/w)*intop_vol(oyyy_e)');
model.component('comp1').variable('var10').set('oyyz_e', 'Jy*(((y*z)/3)-((1/5)*rsqr*delta0))+y*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))+Jy*(((y*z)/3)-((1/5)*rsqr*delta0))+y*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))+Jz*(((y*y)/3)-((1/5)*rsqr*delta1))+z*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Oyyz_e', '(i/w)*intop_vol(oyyz_e)');
model.component('comp1').variable('var10').set('oyzx_e', 'Jy*(((z*x)/3)-((1/5)*rsqr*delta0))+y*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))+Jz*(((y*x)/3)-((1/5)*rsqr*delta0))+z*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))+Jx*(((z*y)/3)-((1/5)*rsqr*delta0))+x*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oyzx_e', '(i/w)*intop_vol(oyzx_e)');
model.component('comp1').variable('var10').set('oyzy_e', 'Jy*(((z*y)/3)-((1/5)*rsqr*delta0))+y*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))+Jz*(((y*y)/3)-((1/5)*rsqr*delta1))+z*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))+Jy*(((z*y)/3)-((1/5)*rsqr*delta0))+y*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oyzy_e', '(i/w)*intop_vol(oyzy_e)');
model.component('comp1').variable('var10').set('oyzz_e', 'Jy*(((z*z)/3)-((1/5)*rsqr*delta1))+y*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))+Jz*(((y*z)/3)-((1/5)*rsqr*delta0))+z*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))+Jz*(((z*y)/3)-((1/5)*rsqr*delta0))+z*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Oyzz_e', '(i/w)*intop_vol(oyzz_e)');
model.component('comp1').variable('var10').set('ozxx_e', 'Jz*(((x*x)/3)-((1/5)*rsqr*delta1))+z*((Jx*x)/3+(x*Jx)/3-((2/5)*rdotj*delta1))+Jx*(((z*x)/3)-((1/5)*rsqr*delta0))+x*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))+Jx*(((x*z)/3)-((1/5)*rsqr*delta0))+x*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Ozxx_e', '(i/w)*intop_vol(ozxx_e)');
model.component('comp1').variable('var10').set('ozxy_e', 'Jz*(((x*y)/3)-((1/5)*rsqr*delta0))+z*((Jx*y)/3+(x*Jy)/3-((2/5)*rdotj*delta0))+Jx*(((z*y)/3)-((1/5)*rsqr*delta0))+x*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))+Jy*(((x*z)/3)-((1/5)*rsqr*delta0))+y*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Ozxy_e', '(i/w)*intop_vol(ozxy_e)');
model.component('comp1').variable('var10').set('ozxz_e', 'Jz*(((x*z)/3)-((1/5)*rsqr*delta0))+z*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))+Jx*(((z*z)/3)-((1/5)*rsqr*delta1))+x*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))+Jz*(((x*z)/3)-((1/5)*rsqr*delta0))+z*((Jx*z)/3+(x*Jz)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Ozxz_e', '(i/w)*intop_vol(ozxz_e)');
model.component('comp1').variable('var10').set('ozyx_e', 'Jz*(((y*x)/3)-((1/5)*rsqr*delta0))+z*((Jy*x)/3+(y*Jx)/3-((2/5)*rdotj*delta0))+Jy*(((z*x)/3)-((1/5)*rsqr*delta0))+y*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))+Jx*(((y*z)/3)-((1/5)*rsqr*delta0))+x*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Ozyx_e', '(i/w)*intop_vol(ozyx_e)');
model.component('comp1').variable('var10').set('ozyy_e', 'Jz*(((y*y)/3)-((1/5)*rsqr*delta1))+z*((Jy*y)/3+(y*Jy)/3-((2/5)*rdotj*delta1))+Jy*(((z*y)/3)-((1/5)*rsqr*delta0))+y*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))+Jy*(((y*z)/3)-((1/5)*rsqr*delta0))+y*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Ozyy_e', '(i/w)*intop_vol(ozyy_e)');
model.component('comp1').variable('var10').set('ozyz_e', 'Jz*(((y*z)/3)-((1/5)*rsqr*delta0))+z*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))+Jy*(((z*z)/3)-((1/5)*rsqr*delta1))+y*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))+Jz*(((y*z)/3)-((1/5)*rsqr*delta0))+z*((Jy*z)/3+(y*Jz)/3-((2/5)*rdotj*delta0))');
model.component('comp1').variable('var10').set('Ozyz_e', '(i/w)*intop_vol(ozyz_e)');
model.component('comp1').variable('var10').set('ozzx_e', 'Jz*(((z*x)/3)-((1/5)*rsqr*delta0))+z*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))+Jz*(((z*x)/3)-((1/5)*rsqr*delta0))+z*((Jz*x)/3+(z*Jx)/3-((2/5)*rdotj*delta0))+Jx*(((z*z)/3)-((1/5)*rsqr*delta1))+x*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Ozzx_e', '(i/w)*intop_vol(ozzx_e)');
model.component('comp1').variable('var10').set('ozzy_e', 'Jz*(((z*y)/3)-((1/5)*rsqr*delta0))+z*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))+Jz*(((z*y)/3)-((1/5)*rsqr*delta0))+z*((Jz*y)/3+(z*Jy)/3-((2/5)*rdotj*delta0))+Jy*(((z*z)/3)-((1/5)*rsqr*delta1))+y*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Ozzy_e', '(i/w)*intop_vol(ozzy_e)');
model.component('comp1').variable('var10').set('ozzz_e', 'Jz*(((z*z)/3)-((1/5)*rsqr*delta1))+z*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))+Jz*(((z*z)/3)-((1/5)*rsqr*delta1))+z*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))+Jz*(((z*z)/3)-((1/5)*rsqr*delta1))+z*((Jz*z)/3+(z*Jz)/3-((2/5)*rdotj*delta1))');
model.component('comp1').variable('var10').set('Ozzz_e', '(i/w)*intop_vol(ozzz_e)');
model.component('comp1').variable('var10').set('OOe', 'sqrt(abs(Oxxx_e)^2+abs(Oxxy_e)^2+abs(Oxxz_e)^2+abs(Oxyx_e)^2+abs(Oxyy_e)^2+abs(Oxyz_e)^2+abs(Oxzx_e)^2+abs(Oxzy_e)^2+abs(Oxzz_e)^2+abs(Oyxx_e)^2+abs(Oyxy_e)^2+abs(Oyxz_e)^2+abs(Oyyx_e)^2+abs(Oyyy_e)^2+abs(Oyyz_e)^2+abs(Oyzx_e)^2+abs(Oyzy_e)^2+abs(Oyzz_e)^2+abs(Ozxx_e)^2+abs(Ozxy_e)^2+abs(Ozxz_e)^2+abs(Ozyx_e)^2+abs(Ozyy_e)^2+abs(Ozyz_e)^2+abs(Ozzx_e)^2+abs(Ozzy_e)^2+abs(Ozzz_e)^2)');
model.component('comp1').variable('var10').label('electric octupole moment');

model.component('comp1').variable.create('var16');
model.component('comp1').variable('var16').set('oxxx_t', '(35*x*x*x*rdotj)-(20*rsqr*(x*x*Jx+x*x*Jx+x*x*Jx))+(delta1*delta1+delta1*delta1+delta1*delta1)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxxx_t', '(1/300)*intop_vol(oxxx_t)');
model.component('comp1').variable('var16').set('oxxy_t', '(35*x*x*y*rdotj)-(20*rsqr*(y*x*Jx+x*x*Jy+x*y*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxxy_t', '(1/300)*intop_vol(oxxy_t)');
model.component('comp1').variable('var16').set('oxxz_t', '(35*x*x*z*rdotj)-(20*rsqr*(z*x*Jx+x*x*Jz+x*z*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxxz_t', '(1/300)*intop_vol(oxxz_t)');
model.component('comp1').variable('var16').set('oxyx_t', '(35*x*y*x*rdotj)-(20*rsqr*(x*x*Jy+y*x*Jx+y*x*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxyx_t', '(1/300)*intop_vol(oxyx_t)');
model.component('comp1').variable('var16').set('oxyy_t', '(35*x*y*y*rdotj)-(20*rsqr*(y*x*Jy+y*x*Jy+y*y*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxyy_t', '(1/300)*intop_vol(oxyy_t)');
model.component('comp1').variable('var16').set('oxyz_t', '(35*x*y*z*rdotj)-(20*rsqr*(z*x*Jy+y*x*Jz+y*z*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxyz_t', '(1/300)*intop_vol(oxyz_t)');
model.component('comp1').variable('var16').set('oxzx_t', '(35*x*z*x*rdotj)-(20*rsqr*(x*x*Jz+z*x*Jx+z*x*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxzx_t', '(1/300)*intop_vol(oxzx_t)');
model.component('comp1').variable('var16').set('oxzy_t', '(35*x*z*y*rdotj)-(20*rsqr*(y*x*Jz+z*x*Jy+z*y*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxzy_t', '(1/300)*intop_vol(oxzy_t)');
model.component('comp1').variable('var16').set('oxzz_t', '(35*x*z*z*rdotj)-(20*rsqr*(z*x*Jz+z*x*Jz+z*z*Jx))+(delta1*delta0+delta1*delta0+delta1*delta0)*(x*rsqr*rdotj+4*(rsqr^2)*Jx)');
model.component('comp1').variable('var16').set('Oxzz_t', '(1/300)*intop_vol(oxzz_t)');
model.component('comp1').variable('var16').set('oyxx_t', '(35*y*x*x*rdotj)-(20*rsqr*(x*y*Jx+x*y*Jx+x*x*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyxx_t', '(1/300)*intop_vol(oyxx_t)');
model.component('comp1').variable('var16').set('oyxy_t', '(35*y*x*y*rdotj)-(20*rsqr*(y*y*Jx+x*y*Jy+x*y*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyxy_t', '(1/300)*intop_vol(oyxy_t)');
model.component('comp1').variable('var16').set('oyxz_t', '(35*y*x*z*rdotj)-(20*rsqr*(z*y*Jx+x*y*Jz+x*z*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyxz_t', '(1/300)*intop_vol(oyxz_t)');
model.component('comp1').variable('var16').set('oyyx_t', '(35*y*y*x*rdotj)-(20*rsqr*(x*y*Jy+y*y*Jx+y*x*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyyx_t', '(1/300)*intop_vol(oyyx_t)');
model.component('comp1').variable('var16').set('oyyy_t', '(35*y*y*y*rdotj)-(20*rsqr*(y*y*Jy+y*y*Jy+y*y*Jy))+(delta1*delta1+delta1*delta1+delta1*delta1)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyyy_t', '(1/300)*intop_vol(oyyy_t)');
model.component('comp1').variable('var16').set('oyyz_t', '(35*y*y*z*rdotj)-(20*rsqr*(z*y*Jy+y*y*Jz+y*z*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyyz_t', '(1/300)*intop_vol(oyyz_t)');
model.component('comp1').variable('var16').set('oyzx_t', '(35*y*z*x*rdotj)-(20*rsqr*(x*y*Jz+z*y*Jx+z*x*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyzx_t', '(1/300)*intop_vol(oyzx_t)');
model.component('comp1').variable('var16').set('oyzy_t', '(35*y*z*y*rdotj)-(20*rsqr*(y*y*Jz+z*y*Jy+z*y*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyzy_t', '(1/300)*intop_vol(oyzy_t)');
model.component('comp1').variable('var16').set('oyzz_t', '(35*y*z*z*rdotj)-(20*rsqr*(z*y*Jz+z*y*Jz+z*z*Jy))+(delta1*delta0+delta1*delta0+delta1*delta0)*(y*rsqr*rdotj+4*(rsqr^2)*Jy)');
model.component('comp1').variable('var16').set('Oyzz_t', '(1/300)*intop_vol(oyzz_t)');
model.component('comp1').variable('var16').set('ozxx_t', '(35*z*x*x*rdotj)-(20*rsqr*(x*z*Jx+x*z*Jx+x*x*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozxx_t', '(1/300)*intop_vol(ozxx_t)');
model.component('comp1').variable('var16').set('ozxy_t', '(35*z*x*z*rdotj)-(20*rsqr*(z*z*Jx+x*z*Jz+x*z*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozxy_t', '(1/300)*intop_vol(ozxy_t)');
model.component('comp1').variable('var16').set('ozxz_t', '(35*z*x*z*rdotj)-(20*rsqr*(z*z*Jx+x*z*Jz+x*z*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozxz_t', '(1/300)*intop_vol(ozxz_t)');
model.component('comp1').variable('var16').set('ozyx_t', '(35*z*y*x*rdotj)-(20*rsqr*(x*z*Jy+y*z*Jx+y*x*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozyx_t', '(1/300)*intop_vol(ozyx_t)');
model.component('comp1').variable('var16').set('ozyy_t', '(35*z*y*y*rdotj)-(20*rsqr*(y*z*Jy+y*z*Jy+y*y*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozyy_t', '(1/300)*intop_vol(ozyy_t)');
model.component('comp1').variable('var16').set('ozyz_t', '(35*z*y*z*rdotj)-(20*rsqr*(z*z*Jy+y*z*Jz+y*z*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozyz_t', '(1/300)*intop_vol(ozyz_t)');
model.component('comp1').variable('var16').set('ozzx_t', '(35*z*z*x*rdotj)-(20*rsqr*(x*z*Jz+z*z*Jx+z*x*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozzx_t', '(1/300)*intop_vol(ozzx_t)');
model.component('comp1').variable('var16').set('ozzy_t', '(35*z*z*y*rdotj)-(20*rsqr*(y*z*Jz+z*z*Jy+z*y*Jz))+(delta1*delta0+delta1*delta0+delta1*delta0)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozzy_t', '(1/300)*intop_vol(ozzy_t)');
model.component('comp1').variable('var16').set('ozzz_t', '(35*z*z*z*rdotj)-(20*rsqr*(z*z*Jz+z*z*Jz+z*z*Jz))+(delta1*delta1+delta1*delta1+delta1*delta1)*(z*rsqr*rdotj+4*(rsqr^2)*Jz)');
model.component('comp1').variable('var16').set('Ozzz_t', '(1/300)*intop_vol(ozzz_t)');
model.component('comp1').variable('var16').set('OOt', 'sqrt(abs(Oxxx_t)^2+abs(Oxxy_t)^2+abs(Oxxz_t)^2+abs(Oxyx_t)^2+abs(Oxyy_t)^2+abs(Oxyz_t)^2+abs(Oxzx_t)^2+abs(Oxzy_t)^2+abs(Oxzz_t)^2+abs(Oyxx_t)^2+abs(Oyxy_t)^2+abs(Oyxz_t)^2+abs(Oyyx_t)^2+abs(Oyyy_t)^2+abs(Oyyz_t)^2+abs(Oyzx_t)^2+abs(Oyzy_t)^2+abs(Oyzz_t)^2+abs(Ozxx_t)^2+abs(Ozxy_t)^2+abs(Ozxz_t)^2+abs(Ozyx_t)^2+abs(Ozyy_t)^2+abs(Ozyz_t)^2+abs(Ozzx_t)^2+abs(Ozzy_t)^2+abs(Ozzz_t)^2)');
model.component('comp1').variable('var16').label('toroidal octupole moment');

model.component('comp1').variable.create('var13');
model.component('comp1').variable('var13').set('dipole_constant', '(1/(12*pi))*(k^4)*mu0_const*(c_const^3)');
model.component('comp1').variable('var13').set('quadrupole_constant', '(1/(160*pi))*(k^6)*mu0_const*(c_const^3)');
model.component('comp1').variable('var13').set('octupole_constant', '(1/(3780*pi))*(k^8)*mu0_const*(c_const^3)');
model.component('comp1').variable('var13').set('mag_dipole_constant', '(1/(12*pi))*(k^4)*mu0_const*(c_const)');
model.component('comp1').variable('var13').set('mag_quadrupole_constant', '(1/(160*pi))*(k^6)*mu0_const*(c_const)');
model.component('comp1').variable('var13').label('constants');

model.component('comp1').variable.create('var14');
model.component('comp1').variable('var14').set('pconj_T', 'conj(px)*((1i*k)/c_const)*Tx + conj(py)*((1i*k)/c_const)*Ty + conj(pz)*((1i*k)/c_const)*Tz');
model.component('comp1').variable('var14').set('Tconj_p', 'conj(((1i*k)/c_const)*Tx)*px + conj(((1i*k)/c_const)*Ty)*py + conj(((1i*k)/c_const)*Tz)*pz');
model.component('comp1').variable('var14').set('pT', 'pconj_T+Tconj_p');

model.component('comp1').variable('var14').set('pconj_T1', 'conj(px)*((1i*(k^3))/c_const)*Tx1 + conj(py)*((1i*(k^3))/c_const)*Ty1 + conj(pz)*((1i*(k^3))/c_const)*Tz1');
model.component('comp1').variable('var14').set('T1conj_p', 'conj(((1i*(k^3))/c_const)*Tx1)*px + conj(((1i*(k^3))/c_const)*Ty1)*py + conj(((1i*(k^3))/c_const)*Tz1)*pz');
model.component('comp1').variable('var14').set('pT1', 'pconj_T1+T1conj_p');

model.component('comp1').variable('var14').set('Tconj_T1', 'conj(((1i*k)/c_const)*Tx)*((1i*(k^3))/c_const)*Tx1 + conj(((1i*k)/c_const)*Ty)*((1i*(k^3))/c_const)*Ty1 + conj(((1i*k)/c_const)*Tz)*((1i*(k^3))/c_const)*Tz1');
model.component('comp1').variable('var14').set('T1conj_T', 'conj(((1i*(k^3))/c_const)*Tx1)*(((1i*k)/c_const)*Tx) + conj(((1i*(k^3))/c_const)*Ty1)*(((1i*k)/c_const)*Ty) + conj(((1i*(k^3))/c_const)*Tz1)*(((1i*k)/c_const)*Tz)');
model.component('comp1').variable('var14').set('TT1', 'Tconj_T1+T1conj_T');

model.component('comp1').variable('var14').set('mconj_Tm', 'conj(mx)*((1i*k)/c_const)*Tx_m + conj(my)*((1i*k)/c_const)*Ty_m + conj(mz)*((1i*k)/c_const)*Tz_m');
model.component('comp1').variable('var14').set('Tmconj_m', 'conj(((1i*k)/c_const)*Tx_m)*mx + conj(((1i*k)/c_const)*Ty_m)*my + conj(((1i*k)/c_const)*Tz_m)*mz');
model.component('comp1').variable('var14').set('mTm', 'mconj_Tm+Tmconj_m');

model.component('comp1').variable('var14').set('Qtconj_Qe', 'conj(Qxx_e)*((1i*k)/c_const)*Qxx_t+conj(Qxy_e)*((1i*k)/c_const)*Qxy_t+conj(Qxz_e)*((1i*k)/c_const)*Qxz_t+conj(Qyx_e)*((1i*k)/c_const)*Qyx_t+conj(Qyy_e)*((1i*k)/c_const)*Qyy_t+conj(Qyz_e)*((1i*k)/c_const)*Qyz_t+conj(Qzx_e)*((1i*k)/c_const)*Qzx_t+conj(Qzy_e)*((1i*k)/c_const)*Qzy_t+conj(Qzz_e)*((1i*k)/c_const)*Qzz_t');
model.component('comp1').variable('var14').set('Qeconj_Qt', 'Qxx_e*conj(((1i*k)/c_const)*Qxx_t)+Qxy_e*conj(((1i*k)/c_const)*Qxy_t)+Qxz_e*conj(((1i*k)/c_const)*Qxz_t)+Qyx_e*conj(((1i*k)/c_const)*Qyx_t)+Qyy_e*conj(((1i*k)/c_const)*Qyy_t)+Qyz_e*conj(((1i*k)/c_const)*Qyz_t)+Qzx_e*conj(((1i*k)/c_const)*Qzx_t)+Qzy_e*conj(((1i*k)/c_const)*Qzy_t)+Qzz_e*conj(((1i*k)/c_const)*Qzz_t)');
model.component('comp1').variable('var14').set('QeQt', 'Qtconj_Qe+Qeconj_Qt');

model.component('comp1').variable('var14').set('Qtmconj_Qm', 'conj(Qxx_m)*((1i*k)/c_const)*Qxx_tm+conj(Qxy_m)*((1i*k)/c_const)*Qxy_tm+conj(Qxz_m)*((1i*k)/c_const)*Qxz_tm+conj(Qyx_m)*((1i*k)/c_const)*Qyx_tm+conj(Qyy_m)*((1i*k)/c_const)*Qyy_tm+conj(Qyz_m)*((1i*k)/c_const)*Qyz_tm+conj(Qzx_m)*((1i*k)/c_const)*Qzx_tm+conj(Qzy_m)*((1i*k)/c_const)*Qzy_tm+conj(Qzz_m)*((1i*k)/c_const)*Qzz_tm');
model.component('comp1').variable('var14').set('Qmconj_Qtm', 'Qxx_m*conj(((1i*k)/c_const)*Qxx_tm)+Qxy_m*conj(((1i*k)/c_const)*Qxy_tm)+Qxz_m*conj(((1i*k)/c_const)*Qxz_tm)+Qyx_m*conj(((1i*k)/c_const)*Qyx_tm)+Qyy_m*conj(((1i*k)/c_const)*Qyy_tm)+Qyz_m*conj(((1i*k)/c_const)*Qyz_tm)+Qzx_m*conj(((1i*k)/c_const)*Qzx_tm)+Qzy_m*conj(((1i*k)/c_const)*Qzy_tm)+Qzz_m*conj(((1i*k)/c_const)*Qzz_tm)');
model.component('comp1').variable('var14').set('QmQtm', 'Qtmconj_Qm+Qmconj_Qtm');

model.component('comp1').variable('var14').set('OTconjOE', 'conj(Oxxx_e)*((1i*k)/c_const)*Oxxx_t+conj(Oxxy_e)*((1i*k)/c_const)*Oxxy_t+conj(Oxxz_e)*((1i*k)/c_const)*Oxxz_t+conj(Oxyx_e)*((1i*k)/c_const)*Oxyx_t+conj(Oxyy_e)*((1i*k)/c_const)*Oxyy_t+conj(Oxyz_e)*((1i*k)/c_const)*Oxyz_t+conj(Oxzx_e)*((1i*k)/c_const)*Oxzx_t+conj(Oxzy_e)*((1i*k)/c_const)*Oxzy_t+conj(Oxzz_e)*((1i*k)/c_const)*Oxzz_t+conj(Oyxx_e)*((1i*k)/c_const)*Oyxx_t+conj(Oyxy_e)*((1i*k)/c_const)*Oyxy_t+conj(Oyxz_e)*((1i*k)/c_const)*Oyxz_t+conj(Oyyx_e)*((1i*k)/c_const)*Oyyx_t+conj(Oyyy_e)*((1i*k)/c_const)*Oyyy_t+conj(Oyyz_e)*((1i*k)/c_const)*Oyyz_t+conj(Oyzx_e)*((1i*k)/c_const)*Oyzx_t+conj(Oyzy_e)*((1i*k)/c_const)*Oyzy_t+conj(Oyzz_e)*((1i*k)/c_const)*Oyzz_t+conj(Ozxx_e)*((1i*k)/c_const)*Ozxx_t+conj(Ozxy_e)*((1i*k)/c_const)*Ozxy_t+conj(Ozxz_e)*((1i*k)/c_const)*Ozxz_t+conj(Ozyx_e)*((1i*k)/c_const)*Ozyx_t+conj(Ozyy_e)*((1i*k)/c_const)*Ozyy_t+conj(Ozyz_e)*((1i*k)/c_const)*Ozyz_t+conj(Ozzx_e)*((1i*k)/c_const)*Ozzx_t+conj(Ozzy_e)*((1i*k)/c_const)*Ozzy_t+conj(Ozzz_e)*((1i*k)/c_const)*Ozzz_t');
model.component('comp1').variable('var14').set('OEconjOT', 'Oxxx_e*conj(((1i*k)/c_const)*Oxxx_t)+Oxxy_e*conj(((1i*k)/c_const)*Oxxy_t)+Oxxz_e*conj(((1i*k)/c_const)*Oxxz_t)+Oxyx_e*conj(((1i*k)/c_const)*Oxyx_t)+Oxyy_e*conj(((1i*k)/c_const)*Oxyy_t)+Oxyz_e*conj(((1i*k)/c_const)*Oxyz_t)+Oxzx_e*conj(((1i*k)/c_const)*Oxzx_t)+Oxzy_e*conj(((1i*k)/c_const)*Oxzy_t)+Oxzz_e*conj(((1i*k)/c_const)*Oxzz_t)+Oyxx_e*conj(((1i*k)/c_const)*Oyxx_t)+Oyxy_e*conj(((1i*k)/c_const)*Oyxy_t)+Oyxz_e*conj(((1i*k)/c_const)*Oyxz_t)+Oyyx_e*conj(((1i*k)/c_const)*Oyyx_t)+Oyyy_e*conj(((1i*k)/c_const)*Oyyy_t)+Oyyz_e*conj(((1i*k)/c_const)*Oyyz_t)+Oyzx_e*conj(((1i*k)/c_const)*Oyzx_t)+Oyzy_e*conj(((1i*k)/c_const)*Oyzy_t)+Oyzz_e*conj(((1i*k)/c_const)*Oyzz_t)+Ozxx_e*conj(((1i*k)/c_const)*Ozxx_t)+Ozxy_e*conj(((1i*k)/c_const)*Ozxy_t)+Ozxz_e*conj(((1i*k)/c_const)*Ozxz_t)+Ozyx_e*conj(((1i*k)/c_const)*Ozyx_t)+Ozyy_e*conj(((1i*k)/c_const)*Ozyy_t)+Ozyz_e*conj(((1i*k)/c_const)*Ozyz_t)+Ozzx_e*conj(((1i*k)/c_const)*Ozzx_t)+Ozzy_e*conj(((1i*k)/c_const)*Ozzy_t)+Ozzz_e*conj(((1i*k)/c_const)*Ozzz_t)');
model.component('comp1').variable('var14').set('OEOT', 'OTconjOE+OEconjOT');
model.component('comp1').variable('var14').label('anapoles');

model.component('comp1').variable.create('var8');
model.component('comp1').variable('var8').set('I_p', 'dipole_constant*(abs(px)^2+abs(py)^2+abs(pz)^2)');

model.component('comp1').variable('var8').set('I_t', 'dipole_constant*(abs(((1i*k)/c_const)*Tx)^2+abs(((1i*k)/c_const)*Ty)^2+abs(((1i*k)/c_const)*Tz)^2)');

model.component('comp1').variable('var8').set('I_T1', 'dipole_constant*(abs(((1i*(k^3))/c_const)*Tx1)^2+abs(((1i*(k^3))/c_const)*Ty1)^2+abs(((1i*(k^3))/c_const)*Tz1)^2)');


model.component('comp1').variable('var8').set('I_pT', 'dipole_constant*(pT)');

model.component('comp1').variable('var8').set('I_TT1', 'dipole_constant*(TT1)');

model.component('comp1').variable('var8').set('I_pT1', 'dipole_constant*(pT1)');

model.component('comp1').variable('var8').set('I_m', 'mag_dipole_constant*(abs(mx)^2+abs(my)^2+abs(mz)^2)');

model.component('comp1').variable('var8').set('I_tm', 'mag_dipole_constant*(abs(((1i*k)/c_const)*Tx_m)^2+abs(((1i*k)/c_const)*Ty_m)^2+abs(((1i*k)/c_const)*Tz_m)^2)');

model.component('comp1').variable('var8').set('I_mTm', 'mag_dipole_constant*(mTm)');

model.component('comp1').variable('var8').set('I_Qe', 'quadrupole_constant*abs(QQe)^2');

model.component('comp1').variable('var8').set('I_Qt', 'quadrupole_constant*abs(((1i*k)/c_const)*QQt)^2');

model.component('comp1').variable('var8').set('I_QeQt', 'quadrupole_constant*(QeQt)');

model.component('comp1').variable('var8').set('I_Qm', 'mag_quadrupole_constant*abs(QQm)^2');

model.component('comp1').variable('var8').set('I_Qtm', 'mag_quadrupole_constant*abs(((1i*k)/c_const)*QQ_tm)^2');

model.component('comp1').variable('var8').set('I_QmQtm', 'mag_quadrupole_constant*(QmQtm)');

model.component('comp1').variable('var8').set('I_Oe', 'octupole_constant*abs(OOe)^2');

model.component('comp1').variable('var8').set('I_Ot', 'octupole_constant*abs(((1i*k)/c_const)*OOt)^2');

model.component('comp1').variable('var8').set('I_OEOT', 'octupole_constant*(OEOT)');

model.component('comp1').variable('var8').set('I_total', 'I_p+I_t+I_T1+ (I_pT)+(I_pT1)+(I_TT1)+I_m+I_tm+(I_mTm)+I_Qe+I_Qt+(I_QeQt)+I_Qm+I_Qtm+(I_QmQtm)+I_Oe+I_Ot+(I_OEOT)');

model.component('comp1').variable('var8').label('scattered power from multipoles');

model.component('comp1').variable.create('var9');
model.component('comp1').variable('var9').set('theta', 'acos(z/r)');
model.component('comp1').variable('var9').set('phi', 'atan2(y,x)');
model.component('comp1').variable('var9').set('r', 'sqrt(x^2+y^2+z^2)');
model.component('comp1').variable('var9').set('Pr', 'intop_surf(ewfd.relPoavx*sin(theta)*cos(phi)+ewfd.relPoavy*sin(theta)*sin(phi)+ewfd.relPoavz*cos(theta))');
model.component('comp1').variable('var9').set('nrelPoav', 'nx*ewfd.relPoavx+ny*ewfd.relPoavy+nz*ewfd.relPoavz');
model.component('comp1').variable('var9').set('scat_pwr', 'intop_nanosurf(nrelPoav)');
model.component('comp1').variable('var9').label('scat power Pr');

model.component('comp1').variable.create('var15');
model.component('comp1').variable('var15').set('U', 'intop_vol(ewfd.Wav)');
model.component('comp1').variable('var15').set('Ue', 'intop_vol(ewfd.Wav)');
model.component('comp1').variable('var15').set('Uh', 'intop_vol(ewfd.Wav)');
model.component('comp1').variable('var15').label('Energy density ');


%% Material

model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat3').selection.set([1 2 3 4 5 7 8 9 10]);
model.component('comp1').material('mat3').label('air');
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', {'na' '0' '0' '0' 'na' '0' '0' '0' 'na'});
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});

model.component('comp1').material.create('mat4', 'Common');
model.component('comp1').material('mat4').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat4').selection.set([6]);
model.component('comp1').material('mat4').label('dielectric');
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('n', {'n' '0' '0' '0' 'n' '0' '0' '0' 'n'});
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});


%% Physics

model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
model.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 2);
model.component('comp1').physics('ewfd').feature('sctr1').selection.named('sel8');
model.component('comp1').physics('ewfd').create('ffd1', 'FarFieldDomain', 3);
model.component('comp1').physics('ewfd').feature('ffd1').selection.set([3 4 5 8 9]);
model.component('comp1').physics('ewfd').feature('ffd1').feature('ffc1').selection.named('sel11');
model.component('comp1').physics('ewfd').create('pmc1', 'PerfectMagneticConductor', 2);
model.component('comp1').physics('ewfd').feature('pmc1').selection.set([1 2 6 7 11 14 20 23 27 30 39]);

model.component('comp1').physics('ewfd').prop('ShapeProperty').set('order_electricfield', 1);
model.component('comp1').physics('ewfd').prop('MeshControl').set('SizeControlParameter', 'Wavelength');
model.component('comp1').physics('ewfd').prop('MeshControl').set('PhysicsControlledMeshMinimumWavelength', 'lda0');
model.component('comp1').physics('ewfd').prop('MeshControl').set('ResolveWaveInLossyMedia', true);
model.component('comp1').physics('ewfd').prop('BackgroundField').set('SolveFor', 'scatteredField');
model.component('comp1').physics('ewfd').prop('BackgroundField').set('Eb', {'exp(-j*ewfd.k0*z)'; '0'; '0'});
model.component('comp1').physics('ewfd').prop('BackgroundField').set('ktMax', '(2*(sqrt(2*log(10))))/ewfd.w0');
model.component('comp1').physics('ewfd').feature('wee1').set('minput_frequency_src', 'userdef');
model.component('comp1').physics('ewfd').feature('wee1').set('minput_frequency', 'root.freq');
% model.component('comp1').physics('ewfd').feature('dcont1').set('pairDisconnect', true);
% model.component('comp1').physics('ewfd').feature('dcont1').label('Continuity');
model.component('comp1').physics('ewfd').feature('sctr1').set('Order', 'SecondOrder');
model.component('comp1').physics('ewfd').feature('ffd1').feature('ffc1').set('SymmetryType0', 'SymmetryInHx');
model.component('comp1').physics('ewfd').feature('ffd1').feature('ffc1').set('SymmetryInPlane1', true);

%% study

model.study.create('std1');
model.study('std1').create('wave', 'Wavelength');

model.sol.create('sol3');
model.sol('sol3').study('std1');
model.sol('sol3').attach('std1');
model.sol('sol3').create('st1', 'StudyStep');
model.sol('sol3').create('v1', 'Variables');
model.sol('sol3').create('s1', 'Stationary');
model.sol('sol3').feature('s1').create('fc1', 'FullyCoupled');

model.sol('sol3').feature('s1').feature.remove('fcDef');

model.study('std1').feature('wave').set('punit', 'm');
model.study('std1').feature('wave').set('plist', 'lda0');

model.sol('sol3').attach('std1');
model.sol('sol3').feature('st1').label('Compile Equations: Wavelength Domain');
model.sol('sol3').feature('v1').label('Dependent Variables 1.1');
model.sol('sol3').feature('v1').set('clistctrl', {'p1'});
model.sol('sol3').feature('v1').set('cname', {'lambda0'});
model.sol('sol3').feature('v1').set('clist', {'lda0'});
model.sol('sol3').feature('s1').label('Stationary Solver 1.1');
model.sol('sol3').feature('s1').feature('dDef').active(true);
model.sol('sol3').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol3').feature('s1').feature('dDef').set('thresh', 0.01);
model.sol('sol3').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol3').feature('s1').feature('aDef').set('complexfun', true);
model.sol('sol3').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol3').runAll;

%% Global Evaluation
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Global Evaluation 1');
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').set('probetag', 'none');


model.result.numerical('gev1').set('looplevelinput', {'manual'});
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').set('expr', {'I_p' 'I_t' 'I_T1' 'I_pT' 'I_pT1' 'I_TT1' 'I_m' 'I_tm' 'I_mTm' 'I_Qe'  ...
'I_Qt' 'I_QeQt' 'I_Qm' 'I_Qtm' 'I_QmQtm' 'I_Oe' 'I_Ot' 'I_OEOT' 'I_total' 'Pr'  ...
'U' 'scat_pwr'});
model.result.numerical('gev1').set('unit', {'W' 'W' 'W' 'W' 'W' 'W' 'W' 'W' 'W' 'W'  ...
'W' 'W' 'W' 'W' 'W' 'W' 'W' 'W' 'W' 'W'  ...
'J' 'W'});
model.result.numerical('gev1').set('descr', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '' '' '' '' '' ''  ...
'' ''});
model.result.numerical('gev1').setResult;


str = mphtable(model,'tbl1');

tbl_data = str.data;



I_p(jlamda)=tbl_data(2);
I_t(jlamda)=tbl_data(3);
I_T1(jlamda)=tbl_data(4);
I_pT(jlamda)=tbl_data(5);
I_TT1(jlamda)=tbl_data(6);
I_pT1(jlamda)=tbl_data(7);

I_m(jlamda)=tbl_data(8);
I_tm(jlamda)=tbl_data(9);
I_mTm(jlamda)=tbl_data(10);

I_Qe(jlamda)=tbl_data(11);
I_Qt(jlamda)=tbl_data(12);
I_QeQt(jlamda)=tbl_data(13);

I_Qm(jlamda)=tbl_data(14);
I_Qtm(jlamda)=tbl_data(15);
I_QmQtm(jlamda)=tbl_data(16);

I_Oe(jlamda)=tbl_data(17);
I_Ot(jlamda)=tbl_data(18);
I_OEOT(jlamda)=tbl_data(19);


jlamda
lamda
n

I_total(jlamda)=tbl_data(20)
Pr(jlamda)=tbl_data(21)

U(jlamda)=tbl_data(22);
scat_pwr(jlamda)=tbl_data(23)


T(jlamda)=toc

%% reset study/solution 
    model.study.remove('std1');
    model.sol.remove('sol1');

end  
% 
% fname=strcat('multi_far_field_LP_torus_half_sym_model_r_0.25q1_R_0.5q1_rsph_Rplusrplus3lda_n_4_6th_20_points_14_07_ldacontmesh.mat');
% save(fname,'n','nnu','llamda','c','q1','small_r','R_1','Pr','I_p','I_m','I_tm','I_mTm','I_t','I_Qe','I_Qm','I_Qtm','I_QmQtm','I_T1','I_Qt','I_pT','I_TT1','I_pT1','I_QeQt','I_Oe','I_Ot','I_OEOT','I_total','U','scat_pwr');

tEnd=toc(tStart)