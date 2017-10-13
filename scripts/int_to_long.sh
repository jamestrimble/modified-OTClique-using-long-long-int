sed -i 's/\bint\b/long long/g' *.c *.h
sed -i 's/long long main/int main/' *.c
sed -i 's/long long argc\b/int argc/' *.c
sed -i 's/%d/%lld/g' *.c
sed -i 's/static long long comp_/static int comp_/' *.c
sed -i 's/INT_MAX/LLONG_MAX/' *.c
