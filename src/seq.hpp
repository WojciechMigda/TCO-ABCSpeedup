/*******************************************************************************
 *
 * Filename: seq.hpp
 *
 * Description:
 *      description
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2015-01-16   xx              Initial version
 *
 ******************************************************************************/

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <utility>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

enum
{
    SEQ_MAXSIZE = 10000
};

bool operator<(std::pair<int, std::string>& p1, std::pair<int, std::string>& p2) {
    if (p1.first != p2.first) {
        return p1.first < p2.first;
    }
    return p1.second.compare(p2.second) < 0;
}

struct Seq {
    std::size_t m_clustno;
    std::string id;
    //std::string v_fam;
    std::string v_gene;
    std::string v_all;
    std::string j_gene;
    std::string j_all;
    std::string junc;
//    std::vector<std::pair<int, std::string>> muts;
    std::vector<std::pair<int, int>> muts;

    int hash_mut(const std::string & str) const
    {
        const char * p = str.c_str();
        return p[0] * 65536 + p[1] * 256 + p[2];
    }

    Seq(std::size_t clustno, jsonxx::Object& o, const std::string& junc_query)
    :
        m_clustno(clustno)
    {
        id = o.get<std::string>("seq_id");
        //v_fam = o.get<jsonxx::Object>("v_gene").get<std::string>("fam");
        v_gene = o.get<jsonxx::Object>("v_gene").get<std::string>("gene");
        v_all = o.get<jsonxx::Object>("v_gene").get<std::string>("all");
        j_gene = o.get<jsonxx::Object>("j_gene").get<std::string>("gene");
        j_all = o.get<jsonxx::Object>("j_gene").get<std::string>("all");
        junc = o.get<std::string>(junc_query);
        if (o.has<jsonxx::Object>("var_muts_nt")) {
            jsonxx::Array a = o.get<jsonxx::Object>("var_muts_nt").get<jsonxx::Array>("muts");
            muts.reserve(a.size());
            for (size_t i = 0; i < a.size(); i++) {
                muts.emplace_back(atoi(a.get<jsonxx::Object>(i).get<std::string>("loc").c_str()),
                    hash_mut(a.get<jsonxx::Object>(i).get<std::string>("mut"))
                );
            }
        }
        std::sort(muts.begin(), muts.end());
    }

//    Seq() {}
};

//bool serialize(Seq& seq, size_t* len, char(* data)[SEQ_MAXSIZE]) {
//    *len = sizeof(size_t) * 7 +
//        seq.id.length() +
//        seq.v_gene.length() +
//        seq.v_all.length() +
//        seq.j_gene.length() +
//        seq.j_all.length() +
//        seq.junc.length();
//    for (size_t i = 0; i < seq.muts.size(); i++) {
//        *len += sizeof(int) + sizeof(size_t) + seq.muts[i].second.length();
//    }
//    if (*len > SEQ_MAXSIZE) {
//        fprintf(stderr, "Seq too long\n");
//        return false;
//    }
//    size_t pos = 0;
//    *(size_t*)((*data) + pos) = seq.id.length();
//    pos += sizeof(size_t);
//    strcpy((*data) + pos, seq.id.c_str());
//    pos += seq.id.length();
//    *(size_t*)((*data) + pos) = seq.v_gene.length();
//    pos += sizeof(size_t);
//    strcpy((*data) + pos, seq.v_gene.c_str());
//    pos += seq.v_gene.length();
//    *(size_t*)((*data) + pos) = seq.v_all.length();
//    pos += sizeof(size_t);
//    strcpy((*data) + pos, seq.v_all.c_str());
//    pos += seq.v_all.length();
//    *(size_t*)((*data) + pos) = seq.j_gene.length();
//    pos += sizeof(size_t);
//    strcpy((*data) + pos, seq.j_gene.c_str());
//    pos += seq.j_gene.length();
//    *(size_t*)((*data) + pos) = seq.j_all.length();
//    pos += sizeof(size_t);
//    strcpy((*data) + pos, seq.j_all.c_str());
//    pos += seq.j_all.length();
//    *(size_t*)((*data) + pos) = seq.junc.length();
//    pos += sizeof(size_t);
//    strcpy((*data) + pos, seq.junc.c_str());
//    pos += seq.junc.length();
//    *(size_t*)((*data) + pos) = seq.muts.size();
//    pos += sizeof(size_t);
//    for (size_t i = 0; i < seq.muts.size(); i++) {
//        *(int*)((*data) + pos) = seq.muts[i].first;
//        pos += sizeof(int);
//        *(size_t*)((*data) + pos) = seq.muts[i].second.length();
//        pos += sizeof(size_t);
//        strcpy((*data) + pos, seq.muts[i].second.c_str());
//        pos += seq.muts[i].second.length();
//    }
//    if (pos != *len) {
//        fprintf(stderr, "serialize error\n");
//        return false;
//    }
//    return true;
//}
//
//void deserialize(char* data, Seq* seq) {
//    size_t pos = 0;
//    size_t len;
//
//    len = *(size_t*)(data + pos);
//    pos += sizeof(size_t);
//    seq->id.clear();
//    for (size_t i = 0; i < len; i++) {
//        seq->id.push_back(*(data + pos++));
//    }
//    len = *(size_t*)(data + pos);
//    pos += sizeof(size_t);
//    seq->v_gene.clear();
//    for (size_t i = 0; i < len; i++) {
//        seq->v_gene.push_back(*(data + pos++));
//    }
//    len = *(size_t*)(data + pos);
//    pos += sizeof(size_t);
//    seq->v_all.clear();
//    for (size_t i = 0; i < len; i++) {
//        seq->v_all.push_back(*(data + pos++));
//    }
//    len = *(size_t*)(data + pos);
//    pos += sizeof(size_t);
//    seq->j_gene.clear();
//    for (size_t i = 0; i < len; i++) {
//        seq->j_gene.push_back(*(data + pos++));
//    }
//    len = *(size_t*)(data + pos);
//    pos += sizeof(size_t);
//    seq->j_all.clear();
//    for (size_t i = 0; i < len; i++) {
//        seq->j_all.push_back(*(data + pos++));
//    }
//    len = *(size_t*)(data + pos);
//    pos += sizeof(size_t);
//    seq->junc.clear();
//    for (size_t i = 0; i < len; i++) {
//        seq->junc.push_back(*(data + pos++));
//    }
//    len = *(size_t*)(data + pos);
//    pos += sizeof(size_t);
//    seq->muts.clear();
//    for (size_t i = 0; i < len; i++) {
//        int loc = *(int *)(data + pos);
//        pos += sizeof(int);
//        std::string mut;
//        size_t len_mut = *(size_t*)(data + pos);
//        pos += sizeof(size_t);
//        for (size_t j = 0; j < len_mut; j++) {
//            mut.push_back(*(data + pos++));
//        }
//        seq->muts.emplace_back(loc, mut);
//    }
//}

#endif /* SEQ_HPP_ */
