#pragma once

namespace visual {
	// calls matplotlib's/plotly's trisurf fuction with data loaded
	// from specified tsv file names
	void trisurf(std::string name_triang, std::string name_sol);

	// calls matplotlib's/plotly's tricotour fuction with data loaded
	// from specified tsv file names
	void tricontour(std::string name_triang, std::string name_sol);
}