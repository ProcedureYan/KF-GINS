#ifndef ODOFILELOADER_H
#define ODOFILELOADER_H

#include "common/types.h"
#include "fileloader.h"

class ODOFileLoader : public FileLoader {

public:
    ODOFileLoader() = delete;
    explicit ODOFileLoader(const string &filename, int columns = 2) {
        open(filename, columns, FileLoader::TEXT);
    }

    const ODO &next() {
        data_ = load(); // 加载新的一条数据
        odo_.time = data_[0];
        odo_.odo_vel = data_[1];
        return odo_;
    }

private:
    ODO odo_;
    vector<double> data_;
};



#endif // ODOFILELOADER_H