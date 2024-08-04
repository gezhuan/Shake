/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <ttp://www.gnu.org/licenses/>  *
 ************************************************************************/

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

struct UserData
{
    double Amp;
    double Tper;
    double dt;
    std::ofstream      oss_ss;       ///< file for stress strain data
};


void Report (DEM::Domain & d2, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (d2.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",d2.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Nc" << Util::_8s << "Nsc \n";
    }
    if (!d2.Finished) 
    {
        size_t Nc = 0;
        size_t Nsc = 0;
        for (size_t i=0; i<d2.CInteractons.Size(); i++)
        {
            Nc += d2.CInteractons[i]->Nc;    //Number of contact
            Nsc += d2.CInteractons[i]->Nsc;  // Number of sliding contact
        }
        dat.oss_ss << Util::_10_6 << d2.Time << Util::_8s << Nc << Util::_8s << Nsc << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

void Setup (DEM::Domain & d2, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //double vel = 0.5*dat.str*std::min(10.0*(dom.Time-dat.T0)/(dat.Tf-dat.T0),1.0);

    double vmod=dat.Amp*M_PI/2.0*sin(d2.Time*M_PI*2.0/dat.Tper);
    d2.GetParticle(-15)->v=Vec3_t(vmod,0.0,0.0);
    d2.GetParticle(-15)->InitializeVelocity(dat.dt);
    d2.GetParticle(-16)->v=Vec3_t(vmod,0.0,0.0);
    d2.GetParticle(-16)->InitializeVelocity(dat.dt);
    d2.GetParticle(-17)->v=Vec3_t(vmod,0.0,0.0);
    d2.GetParticle(-17)->InitializeVelocity(dat.dt);
    d2.GetParticle(-18)->v=Vec3_t(vmod,0.0,0.0);
    d2.GetParticle(-18)->InitializeVelocity(dat.dt);
    d2.GetParticle(-19)->v=Vec3_t(vmod,0.0,0.0);
    d2.GetParticle(-19)->InitializeVelocity(dat.dt);
    d2.GetParticle(-20)->v=Vec3_t(vmod,0.0,0.0);
    d2.GetParticle(-20)->InitializeVelocity(dat.dt);
 
}



int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    
    String ptype;
    bool   Render = true;
    double L_x       ;
    double L_y       ;
    double L_z       ;
    double nu        ;
    double dx        ;
    double dt        ;
    double Dp        ;
    double R         ;   
    double Tf1       ;
    double Tf2       ;
    double dtOut1    ;
    double dtOut2    ;
    double height    ;
    double rho_s     ;
    double rho_l     ;
    double Mu        ;
    double Muw       ;
    double Beta      ;        //
	double Kn        ;         // Normal stiffness
	double Kt        ;         // Tangential stiffness
	double Gn        ;         // Normal dissipative coefficient
	double ratiots   ;     // Coefficient to modify the time step
	double ratiokn   ;
    double R_random  ;
    double Gnw       ;
    double sizedis    ;
    double MUP       ;
    double Ang       ;
    double fraction  ;
    int    Nump      ;
    double Amp       ;
    double Tper      ;



    
    {
        infile >> ptype;           infile.ignore(200,'\n');
        infile >> Render;          infile.ignore(200,'\n');
        infile >> L_x;             infile.ignore(200,'\n');
        infile >> L_y;             infile.ignore(200,'\n');
        infile >> L_z;             infile.ignore(200,'\n');
        infile >> nu;              infile.ignore(200,'\n');
        infile >> dx;              infile.ignore(200,'\n');
        infile >> dt;              infile.ignore(200,'\n');
        infile >> Dp;              infile.ignore(200,'\n');
        infile >> R;               infile.ignore(200,'\n');
        infile >> Tf1;             infile.ignore(200,'\n');
        infile >> Tf2;             infile.ignore(200,'\n');
        infile >> dtOut1;          infile.ignore(200,'\n');
        infile >> dtOut2;          infile.ignore(200,'\n');
        infile >> height;          infile.ignore(200,'\n');
        infile >> rho_s;           infile.ignore(200,'\n');
        infile >> rho_l;           infile.ignore(200,'\n');
        infile >> Mu;              infile.ignore(200,'\n');
        infile >> Muw;             infile.ignore(200,'\n');
        infile >> Beta;            infile.ignore(200,'\n');
        infile >> Kn;              infile.ignore(200,'\n');
        infile >> Gn;              infile.ignore(200,'\n');
        infile >> ratiots;         infile.ignore(200,'\n');
        infile >> ratiokn;         infile.ignore(200,'\n');
        infile >> R_random;        infile.ignore(200,'\n');
        infile >> Gnw;             infile.ignore(200,'\n');
        infile >> sizedis;          infile.ignore(200,'\n');
        infile >> MUP;             infile.ignore(200,'\n');
        infile >> Ang;             infile.ignore(200,'\n');
        infile >> fraction;             infile.ignore(200,'\n');
        infile >> Nump;             infile.ignore(200,'\n');
        infile >> Amp;             infile.ignore(200,'\n');
        infile >> Tper;             infile.ignore(200,'\n');

    }

    Kt=ratiokn*Kn;
  

    double acc=Dp * (rho_s - rho_l) / rho_s;
    double ang=Ang*M_PI/180.0;
    cout <<"acc"<< acc<<endl;
    cout <<"L_z"<< L_z<<endl;
    double Sphere_size;
    double dtdem;
    Vec3_t Xmin0;
    Vec3_t Xmax0;
    Vec3_t Center;
    double R_base=9.2;
    Array<int> delpar0;
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    double thichness=4*R;
    


    DEM::Domain d2;
    
    UserData dat;
    d2.UserData = &dat;
    dat.Amp=Amp;
    dat.Tper=Tper;
    Center = Vec3_t(0.0,0.0,0.0);
    
    d2.AddPlane (-17, Vec3_t(Center(0)-L_x/2.0,Center(1),Center(2)), R, L_z, L_y, 1.0, ang, &axis1);
    d2.AddPlane (-18, Vec3_t(Center(0)+L_x/2.0,Center(1),Center(2)), R, L_z, L_y, 1.0, ang, &axis1);
    d2.AddPlane (-19, Vec3_t(Center(0),Center(1)-L_y/2.0,Center(2)), R, L_x, L_z, 1.0, M_PI/2.0, &axis0);
    d2.AddPlane (-20, Vec3_t(Center(0),Center(1)+L_y/2.0,Center(2)), R, L_x, L_z, 1.0, M_PI/2.0, &axis0);
    d2.AddPlane (-15, Vec3_t(Center(0),Center(1),Center(2)-L_z/2.0), R, L_z, L_y, 1.0, 0.0, &axis1);
    d2.AddPlane (-16, Vec3_t(Center(0),Center(1),Center(2)+L_z/2.0), R, L_z, L_y, 1.0, 0.0, &axis1);
    Xmin0=Vec3_t(Center(0)-L_x/2.0+2.0*R, Center(1)-L_y/2.0 +2.0*R , Center(2)-L_z/2.0+2.0*R);
    Xmax0=Vec3_t(Center(0)+L_x/2.0-2.0*R, Center(1)+L_y/2.0 -2.0*R , Center(2)+L_z/2.0-2.0*R);
    d2.GenSpheresBox (-1, Xmin0, Xmax0, R, rho_s, "HCP",  1234, fraction, sizedis);

    d2.GetParticle(-17)->FixVeloc(); 
    d2.GetParticle(-18)->FixVeloc(); 
    d2.GetParticle(-19)->FixVeloc(); 
    d2.GetParticle(-20)->FixVeloc(); 

    d2.GetParticle(-15)->FixVeloc(); 
    d2.GetParticle(-16)->FixVeloc(); 
    size_t count=0;
    for (size_t np=0;np<d2.Particles.Size();np++)
    {
        
        if (d2.Particles[np]->Tag ==-1)
        {
            count=count+1;
        }
                
        if (d2.Particles[np]->Tag ==-1&&count>Nump)
        {
            d2.Particles[np]->Tag = 10;
        }
    }
    Array<int> delpar1;
    delpar1.Push(10);
    if (delpar1.Size()>0) d2.DelParticles(delpar1);
    for (size_t np=0;np<d2.Particles.Size();np++)
    {
        if (d2.Particles[np]->Tag == -1)
        {
        
            d2.Particles[np]->Ff = d2.Particles[np]->Props.m*Vec3_t(0.0,0.0,acc);
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gn; // restitution coefficient
            d2.Particles[np]->Props.Mu = MUP; // frictional coefficient
        }

        if (d2.Particles[np]->Tag <= -15)
        {
            d2.Particles[np]->v=Vec3_t(Amp,0.0,0.0);
            d2.Particles[np]->Props.Kn = Kn; // normal stiffness
            d2.Particles[np]->Props.Kt = Kt; // trangential stiffness
            d2.Particles[np]->Props.Gn = Gn; // restitution coefficient
            d2.Particles[np]->Props.Mu = Muw; // frictional coefficient
        }
    }

    dtdem = ratiots*d2.CriticalDt(); //Calculating time step
    dat.dt=dtdem;
    d2.Alpha = R/4.0; //Verlet distance
    d2.Solve(/*tf*/Tf2, dtdem, /*dtOut*/dtOut2, Setup, Report, "shakecpu", 2, Nproc);
    d2.Save("Stage_cpu");

}
MECHSYS_CATCH
