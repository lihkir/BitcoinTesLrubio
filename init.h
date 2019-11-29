if (global::idx_test == 1)
{
    std::vector<double> dens;
    dens.push_back(8.590500e-04);
    dens.push_back(6.409600e-03);
    dens.push_back(4.430900e-02);
    dens.push_back(7.928000e-02);
    dens.push_back(4.706500e-02);
    dens.push_back(1.571200e-02);
    dens.push_back(5.720500e-03);
    dens.push_back(1.758300e-03);
    
    for (int i = 0; i < pt_test->M_rows; i++) 
        for (int j = 0 ; j < N_cols; j++)
            u0[i][j]=dens[i];
}
else if (global::idx_test == 2) 
{
    double dx = pt_test->L/N_cols;
    double x;

    for (int i = 0; i < pt_test->M_rows; i++) 
    {
        for (int j = 0 ; j < N_cols; j++)
        {
            x = dx/2 + (i-1)*dx;
            u0[i][j] = 0.12*exp(-200*(x-0.5)*(x-0.5));            
        }
    }  
}
else if ( global::idx_test == 3 ) 
{
    std::vector<double> dens;
    dens.push_back(8.590500e-04);
    dens.push_back(6.409600e-03);
    dens.push_back(4.430900e-02);
    dens.push_back(7.928000e-02);
    dens.push_back(4.706500e-02);
    dens.push_back(1.571200e-02);
    dens.push_back(5.720500e-03);
    dens.push_back(1.758300e-03);
//    dens(1)=3.288742053268562e-03;
//    dens(2)=1.138015958334105e-01;
//    dens(3)=2.500988514694888e-01;
//    dens(4)=9.920555854169366e-02;
//    dens(5)=2.305478842696419e-02;
//    dens(6)=8.206216521787478e-03;
//    dens(7)=5.020573462087361e-03;
//    dens(8)=1.834930402387243e-03;

//    double dx=L/M;
//    double x;

//    for ( int l=1; l<=N; l++ ) {
//      for ( int i=1; i<=M; i++ ) {
//        u0(l, i)=dens(l);
//      }
//    }
}
//  else if ( global::idx_test == 4 ) {

//    grid dens(1, &N);
//    dens(1)=0.2;

//    double dx=L/M;
//    double x;

//    for ( int l=1; l<=N; l++ ) {
//      for ( int i=1; i<=M; i++ ) {
//        u0(l, i)=dens(l);
//      }
//    }
//  }
//  else {
//    std::cerr<<"Test "<<test<<"Not defined"<<std::endl;
//    exit(1);
//  }