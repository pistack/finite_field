template<typename field>
void gauss_elimination<field>::transpose(field* mat, int dim)
{
    int rowi, rowj;
    field tmp;
    int row_idx;
    rowi = 0;
    for(int i=0; i<dim; i++)
    {
        rowj = rowi+dim;
        for(int j=i+1; j<dim; j++)
        {
            tmp = mat[rowi+j];
            mat[rowi+j] = mat[rowj+i]; 
            mat[rowj+i] = tmp;
            rowj += dim;
        }
        rowi += dim;
    }
    return;
}

template<typename field>
field gauss_elimination<field>::reduce(field* mat, int dim)
{
    int rowi, rowj;
    field a_inv, tmp, det(1);
    rowi = 0;
    for(int i=0; i<dim; i++)
    {
        if(mat[rowi+i] == 0)
        {
            rowj = rowi+dim;
            for(int j=i+1; j<dim; j++)
            {
                if(!(mat[rowj+i] == 0))
                {
                    det *= field(-1);
                    for(int k=0; k<dim; k++)
                    {
                        tmp = mat[rowi+k];
                        mat[rowi+k] = mat[rowj+k];
                        mat[rowj+k] = tmp;
                    }
                    break;
                }
                rowj += dim;
            }
        }
        if(mat[rowi+i] == 0)
        return field(0);
        det *= mat[rowi+i];
        a_inv = field(1)/mat[rowi+i];
        for(int k=0; k<dim; k++)
        mat[rowi+k] *= a_inv;
        rowj = rowi+dim;
        for(int j=i+1; j<dim; j++)
        {
            tmp = mat[rowj+i];
            if(!(tmp == 0))
            {
                for(int k=0; k<dim; k++)
                mat[rowj+k] -= tmp*mat[rowi+k];
            }
            rowj += dim;
        }
        rowi += dim;
    }
    return det;
}

template<typename field>
bool gauss_elimination<field>::solve(field* mat, field* y, int dim)
{
    int rowi, rowj;
    field a_inv, tmp;
    rowi = 0;
    for(int i=0; i<dim; i++)
    {
        if(mat[rowi+i] == 0)
        {
            rowj = rowi+dim;
            for(int j=i+1; j<dim; j++)
            {
                if(!(mat[rowj+i] == 0))
                {
                    tmp = y[i]; y[i] = y[j]; y[j] = tmp;
                    for(int k=0; k<dim; k++)
                    {
                        tmp = mat[rowi+k];
                        mat[rowi+k] = mat[rowj+k];
                        mat[rowj+k] = tmp;
                    }
                    break;
                }
                rowj += dim;
            }
        }
        if(mat[rowi+i] == 0)
        return false;
        a_inv = field(1)/mat[rowi+i];
        y[i] *= a_inv;
        for(int k=0; k<dim; k++)
        mat[rowi+k] *= a_inv;
        rowj = rowi+dim;
        for(int j=i+1; j<dim; j++)
        {
            tmp = mat[rowj+i];
            if(!(tmp == 0))
            {
                y[i] -= tmp*y[j];
                for(int k=0; k<dim; k++)
                mat[rowj+k] -= tmp*mat[rowi+k];
            }
            rowj += dim;
        }
        rowi += dim;
    }
    // Solve Mx=y for RREF matrix.
    rowi = dim*(dim-1);
    for(int i=dim-1; i>=0; i--)
    {
        for(int j=dim-1; j>i; j--)
        y[i] -= mat[rowi+j]*y[j];
        rowi -= dim;
    }
    return true;
}

template<typename field>
bool gauss_elimination<field>::inverse(field* mat, field* inv, int dim)
{
    int rowi, rowj;
    field a_inv, tmp1, tmp2;

    for(int i=0; i<dim*dim; i++)
    inv[i] = 0;
    for(int i=0; i<dim*dim; i += dim+1)
    inv[i] = 1;

    rowi = 0;
    for(int i=0; i<dim; i++)
    {
        if(mat[rowi+i] == 0)
        {
            rowj = rowi+dim;
            for(int j=i+1; j<dim; j++)
            {
                if(!(mat[rowj+i] == 0))
                {
                    for(int k=0; k<dim; k++)
                    {
                        tmp1 = mat[rowi+k];
                        mat[rowi+k] = mat[rowj+k];
                        mat[rowj+k] = tmp1;
                        tmp2 = inv[rowi+k];
                        inv[rowi+k] = inv[rowj+k];
                        inv[rowj+k] = tmp2;
                    }
                    break;
                }
                rowj += dim;
            }
        }

        if(mat[rowi+i] == 0)
        return false;

        a_inv = field(1)/mat[rowi+i];
        for(int k=0; k<dim; k++)
        {
            mat[rowi+k] *= a_inv;
            inv[rowi+k] *= a_inv;
        }

        rowj = rowi+dim;
        for(int j=i+1; j<dim; j++)
        {
            tmp1 = mat[rowj+i];
            if(!(tmp1 == 0))
            {
                for(int k=0; k<dim; k++)
                {
                    mat[rowj+k] -= tmp1*mat[rowi+k];
                    inv[rowj+k] -= tmp1*inv[rowi+k];
                }
            }
            rowj += dim;
        }
        rowi += dim;
    }
    // Inverting RREF matrix.
    rowi = dim*(dim-1);
    for(int i=dim-1; i>=0; i--)
    {
        rowj = dim*(dim-1);
        for(int j=dim-1; j>i; j--)
        {
            for(int k=0; k<dim; k++)
            inv[rowi+k] -= mat[rowi+j]*inv[rowj+k];
            rowj -= dim;
        }
        rowi -= dim;
    }
    return true;
}

