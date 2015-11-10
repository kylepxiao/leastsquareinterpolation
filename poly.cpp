#include "poly.h"

poly::poly()
{
    std::vector<double> vec;
    coord.push_back(vec);
    coord.push_back(vec);
}

poly::~poly()
{
    //dtor
}

void poly::getpoints()
{
    coord[0].clear();
    coord[1].clear();
    const char* space = " ";
    bool run = true;
    double num1, num2;
    int pos;
    std::string ans, constans;
    std::string str1, str2;
    cout<<"Enter the degree for which you would like to solve."<<endl;
    cin>>degree;
    cout<<"Enter the points. Put a space in between each value."<<endl;
    std::getline(std::cin, ans);
    while(run)
    {
        std::getline(std::cin, ans);
        constans = ans;
        if (ans.empty())
        {
            run = false;
        }
        else
        {
            for (unsigned int i=0; i<ans.size(); i++)
            {
                if (ans[i] == *space)
                {
                    pos = i;
                    i=ans.size();
                }
            }
            str1 = ans.erase(pos, ans.size()-pos);
            ans = constans;
            str2 = ans.erase(0, (pos+1));
            num1 = atof(str1.c_str());
            num2 = atof(str2.c_str());
            coord[0].push_back(num1);
            coord[1].push_back(num2);
        }
    }
}

std::vector< std::vector<double> > poly::round(std::vector< std::vector<double> > mat)
{
    for (unsigned int i=0; i<mat.size(); i++)
    {
        for (unsigned int k=0; k<mat[i].size(); k++)
        {
            mat[i][k] = floor(mat[i][k]*100+0.5)/100;
        }
    }
    return mat;
}

void poly::printvec(std::vector<double> vec)
{
    for (unsigned int k=0; k<vec.size(); k++)
    {
        cout<<vec[k]<<" ";
    }
    cout<<endl;
}

double poly::dot(std::vector< double > vec1, std::vector< double > vec2)
{
    double num = 0;
    if (!(vec1.size()==vec2.size()))
    {
        cerr<<"Vectors cannot be dot multiplied."<<endl;
        return num;
    }
    for (unsigned int i=0; i<vec1.size(); i++)
    {
        num = num+(vec1[i]*vec2[i]);
    }
    return num;
}

std::vector< std::vector<double> > poly::mult(std::vector< std::vector<double> > mat1, std::vector< std::vector<double> > mat2)
{
    std::vector<double> vec1, vec2;
    std::vector< std::vector<double> > tempmat1;
    if (!(mat1[0].size()==mat2.size()))
    {
        cerr<<"Matrix can not be multiplied"<<endl;
        std::vector< std::vector<double> > vec;
        return vec;
    }
    for (unsigned int i=0; i<mat1.size(); i++)
    {
        vec2.clear();
        for (unsigned int k=0; k<mat2[0].size(); k++)
        {
            vec1.clear();
            for (unsigned int j=0; j<mat2.size(); j++)
            {
                vec1.push_back(mat2[j][k]);
            }
            vec2.push_back(poly::dot(mat1[i], vec1));
        }
        tempmat1.push_back(vec2);
    }
    return tempmat1;
}

double poly::det(std::vector< std::vector<double> > mat)
{
    double tempdet, findet, posval, total=0;
    std::vector< std::vector<double> > tempmat;
    if (!(mat.size()==mat[0].size()))
    {
        cerr<<"Cannot find determinant."<<endl;
        return 0;
    }

    for (unsigned int i=0; i<mat[0].size(); i++)
    {
        posval = mat[0][i];
        tempmat = mat;
        if (mat.size() == 2)
        {
            tempdet = (tempmat[0][0]*tempmat[1][1])-(tempmat[0][1]*tempmat[1][0]);
            return tempdet;
        }
        else
        {
            for (unsigned int k=0; k<mat.size(); k++)
                {
                    tempmat[k].erase(tempmat[k].begin()+i);
                }
            tempmat.erase(tempmat.begin());
            findet = posval*poly::det(tempmat);
            if ((i%2==1))
            {
                findet *= -1;
            }
            total+= findet;
        }
    }
    return total;
}

std::vector< std::vector<double> > poly::invert(std::vector< std::vector<double> > mat)
{
    double vecratio, divideratio;
    int iterate = 0;
    std::vector< std::vector<double> > ident;
    std::vector<double> tempvec, origvec, identvec, origdent, subvector;
    if (!(mat.size()==mat[0].size()))
    {
        cerr<<"Cannot invert matrix."<<endl;
        return ident;
    }
    else if (poly::det(mat)==0)
    {
        cerr<<"Matrix inversion undefined."<<endl;
        return ident;
    }
    for (unsigned int i=0; i<mat.size(); i++)
    {
        tempvec.clear();
        for (unsigned int k=0; k<mat[0].size(); k++)
        {
            if (i==k)
            {
                tempvec.push_back(1);
            }
            else
            {
                tempvec.push_back(0);
            }
        }
        ident.push_back(tempvec);
    }

    for (unsigned int i=0; i<mat.size(); i++)
    {
        tempvec.clear();
        iterate = 1;
        origvec = mat[i];
        origdent = ident[i];
        while (mat[i][i]==0)
        {
            tempvec = mat[i+iterate];
            identvec = ident[i+iterate];
            mat[i] = tempvec;
            ident[i] = identvec;
            mat[i+iterate] = origvec;
            ident[i+iterate] = origdent;
            iterate++;
        }
        divideratio = mat[i][i];
        for (unsigned int k=0; k<mat[0].size(); k++)
        {
            mat[i][k] = mat[i][k]/divideratio;
            ident[i][k] = ident[i][k]/divideratio;
        }
        for (unsigned int j=0; j<mat.size(); j++)
        {
            if (!(j==i))
            {
                tempvec.clear();
                vecratio = mat[j][i];
                for (unsigned int l=0; l<mat[0].size(); l++)
                {
                    mat[j][l] -= vecratio*mat[i][l];
                    ident[j][l] -= vecratio*ident[i][l];
                }
            }
        }
    }
    return ident;
}


//interpolation loop


std::vector< std::vector<double> > poly::interpolate()
{
    double sum;
    std::vector<double> xysums, vec;
    std::vector< std::vector<double> > coeff, constval, tempval, ansmat;
    for (int i=0; i<(degree+1); i++)
    {
        xysums.clear();
        sum = 0;
        for (unsigned int k=0; k<coord[0].size(); k++)
        {
            sum = sum + ((pow(coord[0][k], i)*coord[1][k]));
        }
        xysums.push_back(sum);
        constval.push_back(xysums);
    }
    std::reverse(constval.begin(), constval.end());


    for (int i=0; i<degree+1; i++)//the current variable differentiated
    {
        vec.clear();
        for (int k=0; k<degree+1; k++)//the current term
        {
            sum = 0;
            for (unsigned int j=0; j<coord[0].size(); j++)//the current iteration in the sum
            {
                sum = sum + pow(coord[0][j], ((degree-i)+(degree-k)));//highest-low
            }
            vec.push_back(sum);
        }

        coeff.push_back(vec);
    }

    tempval = poly::invert(coeff);

    if (tempval != ansmat)
    {
        ansmat = poly::mult(poly::invert(coeff), constval);
    }
    else
    {
        cerr<<"The data given is insufficient."<<endl;
    }
    return ansmat;
}

void poly::printans(std::vector< std::vector<double> > mat)
{
    for (unsigned int k=0; k<mat.size(); k++)
    {
        for (unsigned int i=0; i<mat[k].size(); i++)
        {
            cout<<mat[k][i];
            cout<<" ";
        }
        cout<<endl;
    }
}

void poly::inter()
{
    bool running = true;
    string userans, yn;
    cout<<"I will solve for splines given points and a degree.\n"<<endl;
    userans = "inter";
    while (running)
    {
        if (userans == "inter")
        {
            poly::getpoints();
            cout<<"Would you like to round?"<<endl;
            cin>>yn;
            if (yn == "yes")
            {
                cout<<"The coefficients from highest to lowest are:"<<endl;
                poly::printans(poly::round(poly::interpolate()));
            }
            else if (yn == "no")
            {
                cout<<"The coefficients from highest to lowest are:"<<endl;
                poly::printans(poly::interpolate());
            }
            else
            {
                cerr<<"Invalid answer; defaulting to yes."<<endl;
                cout<<"The coefficients from highest to lowest are:"<<endl;
                poly::printans(poly::round(poly::interpolate()));
            }

        }
        else if (userans == "exit")
        {
            running = false;
        }
        else
        {
            cerr<<"Invalid command."<<endl;
        }
        userans = "";
        cout<<"Enter a command to continue."<<endl;
        cin>>userans;
    }
}

