AMREX_FORCE_INLINE
void operator += (const OP_CLASS &rhs)
{
    #define X(name) \
        name += rhs.name;
        OP_VARS
    #undef X
    #ifdef OP_LOOP_VARS
    for (int i = 0; i < OP_LOOP_NUM; i++)
    {
        #define X(name)\
        name[i] += rhs.name[i];
        OP_LOOP_VARS
        #undef X
    }
    #endif
}    

AMREX_FORCE_INLINE
OP_CLASS operator * (const Set::Scalar alpha) const
{
    OP_CLASS ret;
    #define X(name) \
        ret.name = name*alpha;
        OP_VARS
    #undef X
    #ifdef OP_LOOP_VARS
    for (int i = 0; i < OP_LOOP_NUM; i++)
    {
        #define X(name)\
        ret.name[i] = name[i]*alpha;
        OP_LOOP_VARS
        #undef X
    }
    #endif

    return ret;
}    
friend OP_CLASS operator * (const Set::Scalar alpha, const OP_CLASS b);
friend OP_CLASS operator + (const OP_CLASS a, const OP_CLASS b);
friend OP_CLASS operator - (const OP_CLASS a, const OP_CLASS b);
