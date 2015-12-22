#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <stdlib.h>

void makeDir(char name[]){
    struct stat buf;
    int ret;
    char dir[256];
    char mkdir[512];
    
    snprintf(dir,256,"%s",name);
    snprintf(mkdir,512,"mkdir %s",dir);
    
    ret=stat(dir, &buf);
    
    if(ret!=0){
        
        ret=system("ls");
        
        if(ret==0){
            
            ret=system(mkdir);
            
            if(ret==0){
                
                printf("\n\n");
                printf("%sフォルダ作成成功! \n ",dir);
                printf("\n\n ");
                
                ret=system("ls");
                
                if(ret!=0){
                    printf("dirコマンド失敗! \n ");
                }
                
            }else{
                printf("%sフォルダ作成失敗! \n ",dir);
            }
            
        }else{
            printf("dirコマンド失敗! \n ");
        }
    }else{
        printf("%sフォルダが存在します \n ",dir);
    }
}