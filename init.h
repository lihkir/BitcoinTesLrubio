if (global::idx_test == 1)
{
    for (int i = 1; i <= pt_test->M_rows; i++) 
        for (int j = 1 ; j <= N_cols; j++)
            u0(i, j) = 0.04;
}
else
{
    std::cerr << "Test "<< global::idx_test << "Not defined" << std::endl;
    exit(1);
 }