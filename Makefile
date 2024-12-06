.PHONY: all
unit:
	cc -o ./build/tests/unit tests/unit/test_math3d.c -Iinclude

triangle:
	cc -o ./build/tests/triangle tests/triangle/main.c -Iinclude -Iexternals

