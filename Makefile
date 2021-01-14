subdiv:  subdiv.l subdiv.y 
	bison -d subdiv.y
	flex subdiv.l 
	cc -Wall -o $@ subdiv.tab.c lex.yy.c
debug:  subdiv.l subdiv.y 
	bison -v -d subdiv.y --debug
	flex subdiv.l 
	cc -Wall -o subdiv subdiv.tab.c lex.yy.c
