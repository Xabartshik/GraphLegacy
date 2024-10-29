// GraphLegacy.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "GraphLegacy.h"
int main() {
    WeightedGraphLegacy<std::string>::runTests();
    WeightedGraphLegacy<std::string> graph;

    graph.insertVertex("A");
    graph.insertVertex("B");
    graph.insertVertex("C");
    graph.insertVertex("D");

    graph.insertEdge("A", "B", 2.0);
    graph.insertEdge("A", "A", 1.0);
    graph.insertEdge("B", "C", 3.0);
    graph.insertEdge("C", "D", 1.0);
    graph.insertEdge("D", "A", 4.0);

    graph.printVertices();
    graph.printEdges();
    std::cout << "Количество ребер: " << graph.getEdgeCount() << std::endl;
    std::cout << "Количество вершин: " << graph.getVertexCount() << std::endl;
    std::cout << "Минимально вохможный путь между вершинами А и С: " << std::endl;
    graph.findShortestPath("A", "C");
    graph.depthFirstSearch("A");
    graph.breadthFirstSearch("A");
    graph.printAdjacencyMatrix();
    return 0;
}