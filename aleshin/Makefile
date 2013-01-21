aleshin: aleshin.f
	f77 -fast -L/unige/nag/SunOS_5.8_sun4u/lib -o aleshin aleshin.f -lnag
test: aleshin
	(echo 2; echo 5; echo 1 1; echo 1 2; echo 2 1; echo 2 2; echo 3 3) | ./aleshin

