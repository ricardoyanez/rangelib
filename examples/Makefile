
CCFLAGS = -g -std=c99 -Wall

test: clean passage.c rangeair.c
	gcc $(CCFLAGS) passage.c -lrange -lm -o passage
	gcc $(CCFLAGS) rangeair.c -lrange -lm -o rangeair

clean:
	rm -f *~ *.o passage rangeair testRange_C_ACLiC_dict_rdict.pcm testRange_C.*
