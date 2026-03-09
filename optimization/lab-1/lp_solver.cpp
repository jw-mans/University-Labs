#include "lp_solver.h"
#include "matrix.h"
#include "simplex.h"
#include "printer.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

void solveLP()
{
    const double EPS = 1e-7;

    std::cout<<"========================================\n";
    std::cout<<"LINEAR PROGRAMMING SOLVER\n";
    std::cout<<"========================================\n";

    std::vector<std::string> var_names =
    {
        "x1","x2","x3",
        "u4","u5",
        "v4","v5",
        "w1","w2","w3"
    };

    int n_vars = var_names.size();
    int n_constr = 4;

    Matrix A(n_constr,n_vars);

    A(0,0)=2; A(0,1)=1; A(0,2)=3;
    A(1,0)=1; A(1,1)=3; A(1,2)=2;
    A(2,0)=-3;A(2,1)=-2;A(2,2)=-1;
    A(3,0)=2; A(3,1)=1; A(3,2)=2;

    A(0,3)=2;  A(0,4)=1;
    A(1,3)=1;  A(1,4)=2;
    A(2,3)=-3; A(2,4)=-2;
    A(3,3)=1;  A(3,4)=3;

    A(0,5)=-2; A(0,6)=-1;
    A(1,5)=-1; A(1,6)=-2;
    A(2,5)=3;  A(2,6)=2;
    A(3,5)=-1; A(3,6)=-3;

    A(0,7)=1;
    A(1,8)=1;
    A(2,9)=1;

    std::vector<double> b = {14,10,-12,15};

    std::vector<double> c(n_vars,0);

    c[0]=-3;
    c[1]=-2;
    c[2]=-4;

    c[3]=1;
    c[4]=-1;

    c[5]=-1;
    c[6]=1;

    int n_art = n_constr;

    std::vector<std::string> art_names;

    for(int i=0;i<n_art;i++)
        art_names.push_back("y"+std::to_string(i+1));

    int total_vars = n_vars + n_art;

    Matrix tableau(n_constr+1,total_vars+1);

    for(int i=0;i<n_constr;i++)
    {
        for(int j=0;j<n_vars;j++)
            tableau(i,j)=A(i,j);

        tableau(i,total_vars)=b[i];
    }

    for(int i=0;i<n_constr;i++)
        tableau(i,n_vars+i)=1;

    // c_M
    std::vector<double> cM(total_vars,0);
    for(int j=n_vars;j<total_vars;j++)
        cM[j]=1;

    // формируем строку W как в Python
    for(int j=0;j<total_vars;j++)
    {
        double sum=0;

        for(int i=0;i<n_constr;i++)
            sum+=tableau(i,j);

        tableau(n_constr,j)=sum-cM[j];
    }

    // RHS
    double sumB=0;
    for(double bi:b)
        sumB+=bi;

    tableau(n_constr,total_vars)=sumB;

    std::vector<std::string> col_labels;

    for(auto& v:var_names)
        col_labels.push_back(v);

    for(auto& a:art_names)
        col_labels.push_back(a);

    col_labels.push_back("b");

    std::vector<std::string> row_labels = art_names;
    row_labels.push_back("W");

    std::cout<<"Solving Phase I\n";

    SimplexResult phase1 =
        simplexPhase2(tableau,row_labels,col_labels,EPS,"Phase I");

    Matrix tab = phase1.tableau;

    double W = tab(tab.rows-1,tab.cols-1);

    std::cout<<"W = "<<W<<"\n";

    if(std::abs(W)>EPS)
    {
        std::cout<<"Problem infeasible\n";
        return;
    }

    std::cout<<"Removing artificial variables\n";

    int new_cols = n_vars + 1;

    Matrix tab2(n_constr+1,new_cols);

    for(int i=0;i<n_constr;i++)
        for(int j=0;j<n_vars;j++)
            tab2(i,j)=tab(i,j);

    for(int i=0;i<n_constr;i++)
        tab2(i,n_vars)=tab(i,tab.cols-1);

    std::vector<std::string> col2;

    for(auto& v:var_names)
        col2.push_back(v);

    col2.push_back("b");

    std::vector<std::string> row2 = phase1.row_labels;
    row2.back()="F";

    // правильная строка оценок

    std::vector<double> cB;

    for(int i=0;i<n_constr;i++)
    {
        std::string var=row2[i];

        int idx=-1;

        for(int j=0;j<n_vars;j++)
            if(var_names[j]==var)
                idx=j;

        if(idx!=-1)
            cB.push_back(c[idx]);
        else
            cB.push_back(0);
    }

    for(int j=0;j<n_vars;j++)
    {
        double val=0;

        for(int i=0;i<n_constr;i++)
            val+=cB[i]*tab2(i,j);

        tab2(n_constr,j)=val-c[j];
    }

    double rhs=0;

    for(int i=0;i<n_constr;i++)
        rhs+=cB[i]*tab2(i,n_vars);

    tab2(n_constr,n_vars)=rhs;

    std::cout<<"Solving Phase II\n";

    SimplexResult phase2 =
        simplexPhase2(tab2,row2,col2,EPS,"Phase II");

    std::map<std::string,double> sol = phase2.solution;

    std::map<std::string,double> x;

    x["x1"]=sol.count("x1")?sol["x1"]:0;
    x["x2"]=sol.count("x2")?sol["x2"]:0;
    x["x3"]=sol.count("x3")?sol["x3"]:0;

    double u4 = sol.count("u4")?sol["u4"]:0;
    double v4 = sol.count("v4")?sol["v4"]:0;

    double u5 = sol.count("u5")?sol["u5"]:0;
    double v5 = sol.count("v5")?sol["v5"]:0;

    x["x4"]=u4-v4;
    x["x5"]=u5-v5;

    std::cout<<"\nOptimal solution\n";

    for(auto& p:x)
        std::cout<<p.first<<" = "<<p.second<<"\n";

    double F =
        3*x["x1"]+
        2*x["x2"]+
        4*x["x3"]-
        x["x4"]+
        x["x5"];

    std::cout<<"\nF max = "<<F<<"\n";

    // проверка ограничений

    std::cout<<"\nConstraint check\n";

    double c1=2*x["x1"]+x["x2"]+3*x["x3"]+2*x["x4"]+x["x5"];
    double c2=x["x1"]+3*x["x2"]+2*x["x3"]+x["x4"]+2*x["x5"];
    double c3=3*x["x1"]+2*x["x2"]+x["x3"]+3*x["x4"]+2*x["x5"];
    double c4=2*x["x1"]+x["x2"]+2*x["x3"]+x["x4"]+3*x["x5"];

    std::cout<<"c1 = "<<c1<<" <= 14\n";
    std::cout<<"c2 = "<<c2<<" <= 10\n";
    std::cout<<"c3 = "<<c3<<" >= 12\n";
    std::cout<<"c4 = "<<c4<<" = 15\n";

    // двойственная задача

    std::cout<<"\nDual solution\n";

    Matrix finalTab=phase2.tableau;
    std::vector<std::string> cols=phase2.col_labels;

    std::map<std::string,double> y;

    for(int j=0;j<cols.size();j++)
    {
        if(cols[j]=="w1")
            y["y1"]=-finalTab(finalTab.rows-1,j);

        if(cols[j]=="w2")
            y["y2"]=-finalTab(finalTab.rows-1,j);

        if(cols[j]=="w3")
            y["y3"]=-finalTab(finalTab.rows-1,j);
    }

    double y4=(1-y["y1"]-2*y["y2"]+2*y["y3"])/3;

    y["y4"]=y4;

    for(auto&p:y)
        std::cout<<p.first<<" = "<<p.second<<"\n";

    double G=
        14*y["y1"]+
        10*y["y2"]-
        12*y["y3"]+
        15*y["y4"];

    std::cout<<"\nDual value G = "<<G<<"\n";

    std::cout<<"\nDuality check\n";

    std::cout<<"F* = "<<F<<"\n";
    std::cout<<"G* = "<<G<<"\n";
    std::cout<<"|F-G| = "<<std::abs(F-G)<<"\n";

    if(std::abs(F-G)<1e-6)
        std::cout<<"Duality theorem holds\n";
    else
        std::cout<<"Duality theorem violated\n";
}