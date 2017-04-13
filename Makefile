all: app

app: app.cpp midos.cpp
	g++ -std=c++11 app.cpp midos.cpp -o app
	
clean: 
	rm -f ./app
