
#include <cstddef>
#include <stdarg.h>
#include <pthread.h>

struct _thread
{
    pthread_t id;
    pthread_attr_t attrib;
    void **argv;

    _thread() = default;
    _thread(void *(*_task) (void *), std::size_t argc, ...) 
    {
        va_list ap;
        argv = new void*[argc];
        va_start(ap, argc);

        for(std::size_t i = 0; i < argc; ++i){
            argv[i] = va_arg(ap,void*);
        }

        va_end(ap);

        pthread_attr_init(&attrib);
        pthread_create(&id, &attrib, _task, argv);
        
    }
    int join(void **_return)
    {
        return pthread_join(id, _return);
    }

    void cleanUP(void)
    {
        delete[] argv;
    }
};